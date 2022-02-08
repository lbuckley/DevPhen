#VIRTUAL RECIPROCAL TRANSPLANTS

library(cowplot)

#load degree day functions
source("DDFunctions.R")

#model predictions
histogram(la.dat$abs.lat) #lat= 20,35,50 degrees
histogram(la.dat$G) #G=150,350
histogram(la.dat$T0) #T0=

la.dat.sub=la.dat[,c("G","abs.lat","T0","Genus")]

mod1 <- lmer(T0 ~ G*abs.lat+
               (1|Genus), na.action = 'na.omit', REML=FALSE, data = la.dat.sub)

#make new data
G= c(150, 150, 150, 350, 350, 350)
abs.lat= c(20,35,50,20,35,50)
new.dat= as.data.frame(cbind(G, abs.lat))
new.dat$Genus= "Aphis"

new.dat$T0 <- predict(mod1, newdata = new.dat)

#show phenology given fits
#climate change responses

#--------------------------------
#SAVE STATION DATA
#FIND CLOSEST GHCN STATIONS

#Read GGCN database inventory http://www.ncdc.noaa.gov/oa/climate/ghcn-daily/
setwd(paste(fdir,"data/",sep=""))
stations=read.table("ghcnd-inventory_old.txt")
names(stations)=c("Id","Lat","Lon","Var","FirstYear","LastYear")
stations.un= unique(stations$Id)
stats.elev= read.csv("ghcnd-stations.csv")

#Restrict stations to those with recent data and at least 20 years of data
#stations= stations[which(stations$LastYear>2010 & (stations$LastYear-stations$FirstYear)>20 ),]
stations= stations[which(stations$LastYear>2013 & stations$FirstYear<1970),]

#manually drop stations known to have incomplete data
stations= stations[-which(stations$Id %in% c("USC00413060","CA004016640","USC00417622","CA004019040","CA004010984")),]
#CA004019040 pretty good, but look for later dates

#Check have TMIN and TMAX dat
stations= stations[which(stations$Var=="TMAX" | stations$Var=="TMIN"),]
stations= spread(stations,Var, Var)
stations= stations[which(stations$TMAX=="TMAX" & stations$TMIN=="TMIN"),]

stat.coords= cbind(stations$Lon, stations$Lat)

#----------------------------------
#ANALYSIS

gen.dat=new.dat
gen.dat$lon=-105 #choose longitude
gen.dat$lat= gen.dat$abs.lat

j1.all= array(NA, dim=c(1,length(1970:2020),6,6) )
t1.all= array(NA, dim=c(1,length(1970:2020),6,6) )
ngen.all= array(NA, dim=c(1,length(1970:2020),6,6) )

#------
#loop stations
stat.inds= order(gen.dat$lat)[!duplicated(sort(gen.dat$lat))]

for(stat.k in 1:length(stat.inds) ){
  min.dist<- order(spDistsN1(stat.coords, as.numeric(gen.dat[stat.k,c("lon","lat")]), longlat = TRUE))[1:100]
  min.site= stations$Id[min.dist]

  #can't access MX stations, check for US
  inds= grep("MX",min.site)
  if(length(inds)>0) min.site= min.site[-inds]

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
  gen.dat$Id[stat.k]= as.character(min.site[ind])
  gen.dat$st.lat[stat.k]= stations$Lat[min.dist[ind]]
  gen.dat$st.lon[stat.k]= stations$Lon[min.dist[ind]]
  gen.dat$st.elev[stat.k]= stations$elev[min.dist[ind]]

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

  dds= array(NA, dim=c( nrow(gen.dat), length(stat.inds), 20, 2) )

  #CALCULATE DEGREE DAYS
  for(pop.k in 1:nrow(gen.dat) ){

    #dds[,pop.k,stat.k]
    dat$dd= apply( dat[,c("tmin","tmax")], MARGIN=1, FUN=degree.days.mat, LDT=gen.dat$T0[pop.k] )

    dat.dd = dat %>% group_by(year) %>% arrange(j) %>% mutate(cs = cumsum(dd))

    #Estimate phenology across years based on DD accumulations
    #cumsum within groups
    #Egg to adult DD, First date beyond threshold

    js= matrix(NA, 20, length(1970:2020) )
    ts=  matrix(NA, 20, length(1970:2020) )

    for(genk in 1:20){

      phen= dat.dd %>%  group_by(year) %>% slice(which.max(cs> (genk* gen.dat$G[pop.k]) ))

      #replace eroneous values
      phen[which(phen$j==1),"j"]=NA #& phen$cs==0
      #drop years without generation
      phen= phen[!is.na(phen$j),]

      #use only years starting 1970
      years.ind=1970:2020
      year.loop= unique(phen$year)
      year.loop= year.loop[which(year.loop>1969 & year.loop<2021) ]

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

    ngen.all[1,,stat.k,pop.k]= apply(js, MARGIN=2, FUN=function(x)length(which(x>1)) )
    #correct for years without data
    ngen.all[1,which(all.na==TRUE),stat.k,pop.k]=NA

    j1.all[1,,stat.k,pop.k]= js[1,]
    t1.all[1,,stat.k,pop.k]= ts[1,]

  } #end population loop

} #loop station

##SAVE OUTPUT
setwd(paste(fdir,"out_recip/" ,sep=""))
#saveRDS(j1.all, "j1_all.rds")

##READ BACK IN
j1.all= readRDS("j1_all.rds")


#============================================================
#PLOT RECIPROCAL TRANPLANT RESULTS

require(reshape2); require(ggplot2)
library(cowplot)

#Does local adaptation enable more generations, earlier phenology, diff temperature?

#surface plots
#for each genus
#across years, for each station, plot data for all populations

j1.y= j1.all[,,1:6,1:6]
j1.all[1,,stat.k,pop.k]

#melt
j1.m= melt(j1.all, varnames=c("gen","year","station","pop") )
j1.m= na.omit(j1.m)

#add labels
j1.m$lat=c(20,35,50)[j1.m$station]
j1.m$Year= years.ind[j1.m$year]
j1.m$G= years.ind[j1.m$year]
j1.m$T0= round(new.dat$T0,2)[j1.m$pop]
j1.m$G= new.dat$G[j1.m$pop]
j1.m$source.lat= new.dat$abs.lat[j1.m$pop]

#restrict years
j1.y= j1.m#[j1.m$Year %in% c(1970, 2020),]

#plot To and G
new.dat$abs.lat= as.factor(new.dat$abs.lat)
lat.plot= ggplot(data=new.dat, aes(x=G, y = T0, color=abs.lat, group=abs.lat, shape=factor(G)))+
  geom_point(size=2)  +geom_smooth(method=lm, se=FALSE)+
  theme_bw(base_size = 16)+theme(legend.position="none")+
  scale_color_viridis(discrete=TRUE)
#  +labs(color="Absolute latitude (°)")

#plot phenology
j.plot= ggplot(data= j1.y) +
  aes(x = Year, y = value, color=factor(source.lat), lty=factor(lat), shape=factor(G)) +
  geom_line() + geom_point(size=2)+
  ylab("First generation phenology (day of year)")+
  #scale_shape_manual(values=c(21, 22))+
  #scale_fill_gradientn(colours = terrain.colors(5))+
  labs(lty="Latitude (°)",color= "Source latitude (°)", shape="G (°)", fill="T0 (°)")+
  theme_bw(base_size = 16)+ guides(linetype = FALSE)+
  scale_color_viridis(discrete=TRUE, alpha=0.6)+
  geom_text(x=2005, y=163, label="50° Latitude", color="black", size=5)+
  geom_text(x=2005, y=95, label="35° Latitude", color="black", size=5)+
  geom_text(x=2005, y=35, label="20° Latitude", color="black", size=5)

#tile plot
#j.plot= ggplot(data= j1.y) +
#  aes(x = factor(dev), y = factor(lat), z = value, fill = value) +
#  geom_tile() +
#  coord_equal() +
#  scale_fill_distiller(palette="Spectral", na.value="white", name="j")+
#  facet_wrap(~Year)+
#  theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1)) +
#  labs(title=genera[genus.k], x="Source Lat", y="Transplant Lat")

#---------

setwd("/Volumes/GoogleDrive/My Drive/Buckley/work/COISphenology/figures/")
pdf("RecipTran_Dec2021.pdf",height = 6, width = 10)

#plot_grid(lat.plot, j.plot, nrow = 1, rel_widths = c(1.7, 6),labels = c('(a)', '(b)'), label_size = 12)
lat.plot +j.plot +plot_layout(ncol = 2, widths = c(2, 6))+ plot_annotation(tag_levels = 'a')

dev.off()
