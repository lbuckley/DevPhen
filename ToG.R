#see https://github.com/lbuckley/ICBseasonality/blob/master/ICBSeasonality_LocalAdapt.R

fdir= "/Volumes/GoogleDrive/My Drive/Seasonality/"

#Willott and Hassell. 4 species in a location, https://besjournals.onlinelibrary.wiley.com/doi/full/10.1046/j.1365-2435.1998.00180.x

#LOAD libraries
library(ggplot2)
library(dplyr)
library(colorRamps)     # for matlab.like(...)
library(akima) #for interpolations
library(tidyr)
library(sp)
library(rnoaa) #http://recology.info/2015/07/weather-data-with-rnoaa/
library(zoo)
library(fields)
library(patchwork)
library(viridisLite)
library(viridis)
library(lme4)
library(car)
library(sjPlot)
library(broom)

day_of_year<- function(day, format="%Y-%m-%d"){
  day=  as.POSIXlt(day, format=format)
  return(as.numeric(strftime(day, format = "%j")))
}

#LOAD DATA
setwd(paste(fdir,"out/",sep="") )

##READ BACK IN
dddat= read.csv("SeasonalityDatabase_MASTER.csv")

##Restrict to dat with lat / lon
dddat= dddat[which(!is.na(dddat$lon) & !is.na(dddat$lat) ),]

#remove phys outliers
dddat$omit=NA
dddat[which(dddat$BDT.C< (-7)),"omit"]="y" #drops 3
dddat[which(dddat$EADDC>2000),"omit"]="y"  #drops 9

#drop outliers
dddat= dddat[which(is.na(dddat$omit)),]

#group by index
dddat$index= as.factor(dddat$index)

dddat= dddat %>% group_by(index) %>% summarise(Species=head(Species)[1],Order=head(Order)[1],Family=head(Family)[1],Genus=head(Genus)[1],Species.1=head(Species.1)[1],BDT.C=mean(BDT.C), UDT.C=mean(UDT.C),EADDC=mean(EADDC),EEDDC=mean(EEDDC), pest=mean(pest), aquatic=mean(aquatic),pupal=mean(pupal),Location=head(Location)[1],lon=mean(lon),lat=mean(lat), quality=head(quality)[1], parasitoid=head(parasitoid)[1])

dddat= as.data.frame(dddat)

#change trait names
names(dddat)[which(names(dddat)=="BDT.C")]<-"T0"
names(dddat)[which(names(dddat)=="EADDC")]<-"G"

#------------------------
#Count of locations by species

la.spec= dddat %>% group_by(Species) %>% summarise(Npop= length(lon)  )
la.spec= subset(la.spec, la.spec$Npop>4)

la.dat.sp= subset(dddat, dddat$Species %in% la.spec$Species)

#---
#map
library(lattice)
library(rasterVis)
library(maps)
library(mapdata)
library(maptools)

#load worldmap
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf()+geom_point(data=la.dat,aes(x=lon, y = lat, color=Species))+ theme(legend.position = "none")
#---

#get elevation
#library(rgbif)
#elevation= elevation(latitude = la.dat$lat, longitude = la.dat$lon,
#                 elevation_model = "gtopo30",
#                 username ="lbuckley")

#la.dat$elev=elevation$elevation_geonames
#la.dat$elev[la.dat$elev<0]=NA

# write out
#setwd(paste(fdir,"out/",sep="") )
#write.csv(la.dat,"SpeciesLA.csv")

#read in with elevation
setwd(paste(fdir,"out/",sep="") )
la.dat.e= read.csv("SpeciesLA.csv")

#------------------------
# BY GENUS

la.gen= dddat %>% group_by(Genus) %>% summarise(Npop= length(lon)  )
la.gen= subset(la.gen, la.gen$Npop>9) #change?

la.dat= subset(dddat, dddat$Genus %in% la.gen$Genus)

#get elevation
#library(rgbif)
#elevation= elevation(latitude = la.dat$lat, longitude = la.dat$lon,
#                 elevation_model = "gtopo30",
#                 username ="lbuckley")

#la.dat$elev=elevation$elevation_geonames
#la.dat$elev[la.dat$elev<0]=NA

# write out
#setwd(paste(fdir,"out/",sep="") )
#write.csv(la.dat,"GenusLA.csv")

#read in with elevation
#setwd(paste(fdir,"out/",sep="") )
#la.dat= read.csv("GenusLA.csv")

#------
#analysis

#la.dat=la.dat.e #include elevation

#by order
orders= unique(la.dat$Order)[1:4]
#Drop "Thysanoptera" due to limited data
la.dat= la.dat[la.dat$Order %in% orders,]

#stats by lat elev
la.dat$abs.lat= abs(la.dat$lat)

#check distributions
ggplot(la.dat, aes(x=abs.lat)) +
  geom_histogram()+facet_grid(.~Order)

#scale
#la.dat2 <- transform(la.dat,
#                    G_cs=scale(G),
#                    abs.lat_cs=scale(abs.lat))

#models
mod1 <- lmer(T0 ~ G*abs.lat+
               (1|Genus), na.action = 'na.omit', REML=FALSE, data = la.dat)

#mod1 <- lmer(T0 ~ G_cs*abs.lat_cs+
#               (1|Genus), na.action = 'na.omit', REML=FALSE, data = la.dat2)

#mod1 <- lmer(T0 ~ G_cs*log(elev)+
#               (1|Genus), na.action = 'na.omit', REML=FALSE, data = la.dat2)

#by order
ord.k=4
mod1 <- lmer(T0 ~ G_cs*abs.lat_cs+
               (1|Genus), na.action = 'na.omit', REML=FALSE,
             data = la.dat2[la.dat2$Order==orders[ord.k],])

Anova(mod1, type=3)
#summary(mod1)

fixef(mod1)

#get standardized coefficients, see https://stackoverflow.com/questions/25142901/standardized-coefficients-for-lmer-model
stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}

stdCoef.merMod(mod1)

p2 <- plot_model(mod1, type="pred", terms=c("G","abs.lat"), show.data=TRUE)
#p2 <- plot_model(mod1, type="pred", terms=c("G_cs","abs.lat_cs"), show.data=TRUE)
p2

#extract coefficients
coef_st = tidy(mod1,
               effects = "fixed",
               conf.int = TRUE,
               conf.method = "profile")

#significant for orders: Coleoptera, Homoptera, ord.k=1,3
#inverse relationship with T0 and G, steeper at higher latitudes

#----
#plot

gen.lab= la.dat[duplicated(la.dat$Genus)==FALSE,]

#by latitude
p<- ggplot(data=la.dat, aes(x=G, y = T0, color=Genus))+facet_grid(.~Order) +
  theme_bw()+ geom_point(aes(size=abs(lat)))  +geom_smooth(method=lm, se=FALSE)+scale_size(range = c(1, 4))+
  theme(legend.position="bottom")+
  labs(size="Absolute latitude (°)")
#geom_text(data=gen.lab, aes(label=Genus), hjust= -0.2)

#plot out
setwd(paste(fdir,"out/",sep="") )
pdf("Fig1_ToG.pdf", height = 6, width = 10)
p
#T0.plot +G.plot+ plot_layout(ncol=1)
dev.off()

pdf("FigS1_ToG.pdf", height = 8, width = 8)
p2
dev.off()

#doesn't align
#Agrotis: includes cutworm pests
#Ceratitis: pests, can tolerate cool conditions
#shallow: Brevicoryne: cabbage aphids




