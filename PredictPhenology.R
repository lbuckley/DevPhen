#VIRTUAL RECIPROCAL TRANSPLANTS

# read data in
fdir= "/Volumes/GoogleDrive/My Drive/Seasonality/"
setwd(paste(fdir,"out/",sep="") )
la.dat=read.csv("SpeciesLA.csv")

#--------------------------------
#FIND CLOSEST GHCN STATIONS

#Read GGCN database inventory http://www.ncdc.noaa.gov/oa/climate/ghcn-daily/
setwd(paste(fdir,"data/",sep=""))
stations=read.table("ghcnd-inventory.txt")
names(stations)=c("Id","Lat","Lon","Var","FirstYear","LastYear")
stations.un= unique(stations$Id)
stats.elev= read.csv("ghcnd-stations.csv")

#Restrict stations to those with recent data and at least 20 years of data
stations= stations[which(stations$LastYear>2010 & (stations$LastYear-stations$FirstYear)>20 ),] 

#Check have TMIN and TMAX dat
stations= stations[which(stations$Var=="TMAX" | stations$Var=="TMIN"),]
stations= spread(stations,Var, Var)
stations= stations[which(stations$TMAX=="TMAX" & stations$TMIN=="TMIN"),]

stat.coords= cbind(stations$Lon, stations$Lat)

#--------------------------------
#SET UP DATA STORAGE
# phen.dat= array(NA, dim=c(nrow(la.dat), length(1970:2015), 20, 10), dimnames=list(NULL,as.character(1970:2015),as.character(1:20),c("phen","Tmean.e","Tsd.e","T10q.e","Tmean","Tsd","T10q","DaysGen", "Tmean.e.fixed", "Tmean.fixed")) )
# ngens= array(NA, dim=c(nrow(la.dat), length(1970:2015)), dimnames=list(NULL,as.character(1970:2015)))
# phen.fixed= array(NA, dim=c(nrow(la.dat), 20, 8), dimnames=list(NULL,as.character(1:20),c("phen","Tmean.e","Tsd.e","T10q.e","Tmean","Tsd","T10q","DaysGen")) )

#----------------------------------
#ANALYSIS

specs= unique(la.dat$Species)

j1.all= array(NA, dim=c(length(specs),length(2015:2020),100,100) )
t1.all= array(NA, dim=c(length(specs),length(2015:2020),100,100) )
ngen.all= array(NA, dim=c(length(specs),length(2015:2020),100,100) )

for(spec.k in 1:length(specs)){
  
  spec.dat= subset(la.dat, la.dat$Species==specs[spec.k])
  stat.inds= order(spec.dat$lat)[!duplicated(sort(spec.dat$lat))]
  
  for(stat.k in 1:length(stat.inds) ){  
    min.dist<- order(spDistsN1(stat.coords, as.numeric(spec.dat[stat.k,c("lon","lat")]), longlat = TRUE))[1:100]
    min.site= stations$Id[min.dist]
    
    ind=0
    years=NA
    
    while(length(years)<10 & ind<100){
      ind= ind + 1
      tmax=try( ghcnd_search(min.site[ind], var = "TMAX") )
      tmin=try( ghcnd_search(min.site[ind], var = "TMIN") )
      while( (is.null(nrow(tmax$tmax)) | is.null(nrow(tmin$tmin))) & ind<100){ ind=ind+1
      tmax=try( ghcnd_search(min.site[ind], var = "TMAX") )
      tmin=try( ghcnd_search(min.site[ind], var = "TMIN") )
      }
      
      if(!is.null(nrow(tmax$tmax)>0) ){ #CHECK DATA
        #combine tmin and tmax
        match1= match(tmax$tmax$date, tmin$tmin$date)
        is.matched= !is.na(match1)
        
        dat= tmax$tmax
        dat$tmin= NA
        dat$tmin[is.matched]= tmin$tmin$tmin[match1[is.matched]]
        
        #split date 
        date= as.Date(dat$date, "%Y-%m-$d")
        dat$year=as.numeric(format(date, "%Y"))
        dat$month=as.numeric(format(date, "%m"))
        dat$day=as.numeric(format(date, "%d"))
        dat$j= unlist(lapply(date, FUN="day_of_year"))
        
        #Format data
        dat$tmax= as.numeric(dat$tmax)
        dat$tmin= as.numeric(dat$tmin)
        
        dat$tmax[which(dat$tmax==-9999)]= NA
        dat$tmax= dat$tmax/10 #correct for tenths of degrees or mm
        dat$tmin[which(dat$tmin==-9999)]= NA
        dat$tmin= dat$tmin/10 #correct for tenths of degrees or mm
        
        #Catch other NA values
        dat$tmax[which(dat$tmax>200)]= NA
        dat$tmin[which(dat$tmin>200)]= NA
        dat$tmax[which(dat$tmax< -200)]= NA
        dat$tmin[which(dat$tmin< -200)]= NA
        
        ## FIND YEARS WITH NEARLY COMPLETE DATA
        dat.agg= aggregate(dat, list(dat$year),FUN=function(x)length(na.omit(x))  )  
        years= dat.agg$Group.1[which(dat.agg$tmax>300)]
        dat= dat[which(dat$year %in% years),]
        
      } #END CHECK NULL
    } #END WHILE YEARS
    
    #RECORD SITE DATA
    spec.dat$Id[stat.k]= as.character(min.site[ind])
    spec.dat$st.lat[stat.k]= stations$Lat[min.dist[ind]]
    spec.dat$st.lon[stat.k]= stations$Lon[min.dist[ind]]
    spec.dat$st.elev[stat.k]= stations$elev[min.dist[ind]]
    
    #--------------------------------------
    #INTERPOLATE MISSING DATA
    dat$tmin= na.approx(dat$tmin, maxgap=7, na.rm = FALSE)
    dat$tmax= na.approx(dat$tmax, maxgap=7, na.rm = FALSE)
    #Tmean
    dat$tmean= (dat$tmin+dat$tmax)/2
    
    #cut missing data
    dat= dat[which(!is.na(dat$tmin) & !is.na(dat$tmax)),]
    
    #---------------------------
    #TRANSPLANT
    
    dds= array(NA, dim=c( nrow(spec.dat), length(stat.inds), 20, 2) )
    
    #normalize S hemisphere to N hemisphere 
    if(spec.dat[stat.k,"lat"]<0) 
    {inds.jd= which(dat$j>181)
    inds.jj= which(dat$j<182)
    
    dat[inds.jd, "j"] = dat$j[inds.jd] -181
    dat[inds.jj, "j"] = dat$j[inds.jj] +184
    } #end fix s hemi
    
    #CALCULATE DEGREE DAYS
    for(pop.k in 1:nrow(spec.dat) ){
      
      #dds[,pop.k,stat.k]
      dat$dd= apply( dat[,c("tmin","tmax")], MARGIN=1, FUN=degree.days.mat, LDT=spec.dat$T0[pop.k] )
      
      dat.dd = dat %>% group_by(year) %>% arrange(j) %>% mutate(cs = cumsum(dd))
      
      #Estimate phenology across years based on DD accumulations
      #cumsum within groups
      #Egg to adult DD, First date beyond threshold
      
      js= matrix(NA, 20, length(1970:2015) )
      ts=  matrix(NA, 20, length(1970:2015) )
      
      for(speck in 1:20){
        
        phen= dat.dd %>%  group_by(year) %>% slice(which.max(cs> (speck* spec.dat$G[pop.k]) ))
        
        #replace eroneous values
        phen[which(phen$j==1),"j"]=NA #& phen$cs==0
        #drop years without generation
        phen= phen[!is.na(phen$j),]
        
        #use only years starting 1970
        years.ind=1970:2015
        year.loop= unique(phen$year)
        year.loop= year.loop[which(year.loop>1969 & year.loop<2016) ]
        
        if(length(year.loop)>0) for(yeark in 1: length(year.loop)){
          
          year1= year.loop[yeark]
          j= as.numeric(phen[phen$year==year1,"j"])
          
          dat.yr= dat.dd[dat.dd$year==year1,]
          
          #FIND GEN TIME
          j.gs= as.numeric(ifelse(genk==1, dat.yr[which.min(dat.yr$dd>0),"j"], js[genk-1,yeark] ))
          
          #TEMPS: ACROSS GENERATIONS
          ts[genk,yeark]= mean(as.numeric(unlist(dat.yr[dat.yr$j %in% j.gs:j,"tmean"])))
          
          js[genk,yeark]= j
          
        } #end year loop
      } #end GEN LOOP
      
      #Number generations by year
      
      ##DROP?
      #years.match= match(year.loop, years.ind)
      
      all.na= apply(js, MARGIN=1, function(x)all(is.na(x)) )
      
      ngen.all[genus.k,,stat.k,pop.k]= apply(js, MARGIN=2, FUN=function(x)length(which(x>1)) )
      #correct for years without data
      ngen.all[genus.k,which(all.na==TRUE),stat.k,pop.k]=NA
      
      j1.all[genus.k,,stat.k,pop.k]= js[1,]
      t1.all[genus.k,,stat.k,pop.k]= ts[1,]
      
    } #end population loop
    
  } #loop station
} #loop genus

##SAVE OUTPUT
setwd(paste(fdir,"out/" ,sep=""))
saveRDS(phen.dat, "phendat_gen.rds")
saveRDS(phen.fixed, "phenfix_gen.rds")
saveRDS(ngens, "ngens_gen.rds")
saveRDS(dddat, "dddat_media_gen.rds")

##READ BACK IN
# setwd(paste(fdir,"out/",sep="") )
# phen.dat= readRDS("phendat_gen.rds")
# phen.fixed= readRDS("phenfix_gen.rds")
# ngens= readRDS("ngens_gen.rds")
# dddat= readRDS("dddat_media_gen.rds")

#drop omit data
dropi= which( dddat$omit=="y" )
phen.dat= phen.dat[-dropi,,,]
phen.fixed= phen.fixed[-dropi,,]
ngens= ngens[-dropi,]
dddat= dddat[-dropi,]

#============================================================
#PLOT RECIPROCAL TRANPLANT RESULTS

require(reshape2); require(ggplot2)
library(cowplot)

#Does local adaptation enable more generations, earlier phenology, diff temperature?

#surface plots
#for each genus
#across years, for each station, plot data for all populations

setwd("/Volumes/GoogleDrive/My Drive/Buckley/work/ICBSeasonality/figures/LocalAdaptation/") 
pdf("RecipTran.pdf",height = 5, width = 10)

for(genus.k in 1:length(genera)){
  
  gen.dat= subset(la.dat, la.dat$Genus==genera[genus.k])
  lat.ord= order(abs(gen.dat$lat))
  #order lats
  lats= round(gen.dat$lat[lat.ord],2)
  
  ngen.g= ngen.all[genus.k,,lat.ord,lat.ord]
  j1.g= j1.all[genus.k,,lat.ord,lat.ord]
  t1.g= t1.all[genus.k,,lat.ord,lat.ord]
  
  #melt
  ngen.m= melt(ngen.g, varnames=c("year","station","pop") )
  j1.m= melt(j1.g, varnames=c("year","station","pop") )
  t1.m= melt(t1.g, varnames=c("year","station","pop") )
  
  #average across years
  ngen.y= aggregate(ngen.m, by= list(ngen.m$station, ngen.m$pop), FUN=mean, na.rm=TRUE )
  j1.y= aggregate(j1.m, by= list(j1.m$station, j1.m$pop), FUN=mean, na.rm=TRUE )
  t1.y= aggregate(t1.m, by= list(t1.m$station, t1.m$pop), FUN=mean, na.rm=TRUE )
  
  #add lat
  ngen.y$source.lat= abs(lats[ngen.y$pop]);  ngen.y$trans.lat= abs(lats[ngen.y$station])
  j1.y$source.lat= abs(lats[j1.y$pop]);  j1.y$trans.lat= abs(lats[j1.y$station])
  t1.y$source.lat= abs(lats[t1.y$pop]);  t1.y$trans.lat= abs(lats[t1.y$station])
  
  #surface plots # x = source.lat, y = trans.lat
  #ngen       
  n.plot= ggplot(ngen.y) + 
    aes(x = factor(source.lat), y = factor(trans.lat), z = value, fill = value) + 
    geom_tile() + 
    coord_equal() + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="number\ngenerations")+
    theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title=genera[genus.k], x="Source Lat", y="Transplant Lat")
  #, breaks=c(1,2, 5,10,20)) 
  #  +theme_bw(base_size=18)+xlab("T0 (°C)")+ylab("G")+ggtitle(lats[i])+ theme(legend.position="right")+ coord_fixed(ratio = 0.01) 
  
  #j1
  j.plot= ggplot(j1.y) + 
    aes(x = factor(source.lat), y = factor(trans.lat), z = value, fill = value) + 
    geom_tile() + 
    coord_equal() + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="j")+
    theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title=genera[genus.k], x="Source Lat", y="Transplant Lat")
  
  #t1
  t.plot= ggplot(t1.y) + 
    aes(x = factor(source.lat), y = factor(trans.lat), z = value, fill = value) + 
    geom_tile() + 
    coord_equal() + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="temperature") + 
    theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title=genera[genus.k], x="Source Lat", y="Transplant Lat")
  
  p1= plot_grid(n.plot, j.plot, t.plot, labels = c("A", "B","C"), ncol=3)
  print(p1)
  
} #end loop genera

dev.off()