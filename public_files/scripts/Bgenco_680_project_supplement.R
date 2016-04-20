#############################   HEADER    ######################################
#         Brandon Genco Script for 680 paper    2015-12-07                     #
#                                                                              #                
#               Update from script for  presentation: 2015-11-20               #
#                                                                              #
#                 requires two external datasets:                              #
#                                                                              #  
# (1) Smith and Sandwell 30-Plus v.11  30 arc second global topography         #
#     http://topex.ucsd.edu/WWW_html/srtm30_plus.html                          #
#     This require conversion 0-360 longitude -> -180 - 180 longitude:         #
#              & subsettinfg for region of inetrest                            #
#     which can be accomplished using command line GMT tools:                  #
#                                                                              #
#          $ GMT grdedit strmthirty.grd −R-180/180/-90/90 −S                   #
#          $ GMT grdreformat strmthirty.grd -R-75/-60/31.5/42.5 ne.grd         #
#                                                                              #                                                       
# (2) CeDMAR biological database "DatabaseSyndeepCedamar.csv" avaible here:    #
#             http://www2.hawaii.edu/~bmgenco/public_files/                    #
#                                                                              #
#                                                                              #
###           MINMIMIZE SECTIONS (RSTUDIO) TO MAKE SCRIPT EASIER TO READ     ###
#############################    Packages 	 ###################################

rm(list=ls())
library(lattice)
library(marmap)
library(sp)
library(geoR)
library(rgdal)
library(rgeos)
library(raster)
library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)

#library(colorspace) 
#pal<-choose_palette()
# for creating a custom color ramps
# sequential color ramps, which tonally should be 
# accesible to more common types of colorblindness.
#############################   Paths & Working Directories   ###################################

wwd<-(path="/home/brandon/athena/rworkd/Link to temp_R_drop/680/final_project/")
# mac directories
hwd<-(path="/Users/Brandon/Dropbox/temp_R_drop/680/final_project/")

plotsave<-(path="/home/brandon/athena/rworkd/680/figures/")

setwd(wwd)
#############################   Conversion of netcdf to r object   ####################
      #                     Coordinates New England/ NW Atlantic                     #
      #                           -75/-60/31.5/42.5                                  #

# nc<-raster("/home/brandon/athena/ncworkd/ne.grd", varname="z")
# nenc<-as.bathy(nc)
# saveRDS(nenc, file="/home/brandon/Dropbox/temp_R_drop/680/final_project/nenc.rda")
                          

                      ### pick NE to use: above coordinates ###
NE<-readRDS("nenc.rda")

#############################   Points site Data   ###################################
        # athena
#o<-read.csv("/home/brandon/athena/rworkd/Link to temp_R_drop/DEEP/data/DatabaseSyndeepCedamar.csv")
#online:
o<-read.csv("http://www2.hawaii.edu/~bmgenco/public_files/DatabaseSyndeepCedamar.csv")

osites<-unique(o[c("LON_DEG", "LAT_DEG")])
oloc<-cbind(osites$LON_DEG, osites$LAT_DEG)
opts<-SpatialPoints(oloc)
proj4string(opts)<-CRS("+proj=longlat +datum=WGS84")

NEbound<-as.SpatialGridDataFrame(NE)
neopts<-opts[NEbound,]

#############################   Data Setsubsetting NW Atlantic ###################################       


#o<-na.omit(o)
#o$MAX_DEPTH<-as.numeric(paste(o$MAX_DEPTH)) must be down after subsetting
nw<-o[o$LON_DEG> -75 & o$LON_DEG< -60 & o$LAT_DEG> 31.5 & o$LAT_DEG< 42.5,]
#easier to type
h<-nw
#Date conversion!!
sites<-unique(h$EnvSiteID)
h$START_DATE<-as.Date(h$START_DATE, "%d/%m/%Y")

#############################   Latlon->UTM; utm zon 19n    #############################     

hloc<-cbind(h$LON_DEG, h$LAT_DEG)
hpts<-SpatialPoints(hloc)
proj4string(hpts)<-CRS("+proj=longlat +datum=WGS84")
proj4string(hpts)


res <- spTransform(hpts, CRS("+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
h.utm<-as.data.frame(res)

h$LON_DEG<-h.utm$coords.x1
h$LAT_DEG<-h.utm$coords.x2


###     test if it works        
#tsites<-unique(h[c("LON_DEG", "LAT_DEG")])
#tloc<-cbind(tsites$LON_DEG, tsites$LAT_DEG)
#tpts<-SpatialPoints(tloc)
#proj4string(tpts)<-CRS("+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")
#par(mfrow=c(2,2))
#plot(tpts)
#plot(neopts)

#############################   Alpha diversity:    #############################

id<-unique(h$EnvSiteID) 
a<-(1:n_distinct(h$EnvSiteID))
for (i in  1:n_distinct(h$EnvSiteID)){
  a[i]<-sum(h$EnvSiteID==id[i])}
  b<-cbind(id,a)


rm(a)
b<-as.data.frame(b)

names(b)<-c("EnvSiteID", "Diversity")


        # depth & coords & date

h$MAX_DEPTH<-as.numeric(paste(h$MAX_DEPTH))

#### !duplicated command could replace other code in this script 
d<-h[!duplicated(h$EnvSiteID), c("EnvSiteID", "LON_DEG", "LAT_DEG", "START_DATE", "MAX_DEPTH")]
names(d)<-c("EnvSiteID",  "Longitude", "Latitude", "Date" , "Depth")

## site 1652 had date error, listed as 04/31/1966, there are only 30 days in april.. produces NAs 

# diversity site + depth, coords, date:
h.div<-merge(b, d, by="EnvSiteID")
rm(b)
h.div$EnvSiteID<-(1:24)
names(h.div)<-c("Site", "Diversity", "Longitude", "Latitude", "Date" , "Depth")

#############################   Abundance:    #############################
        ## species selection (replace this with random & set seed):
#sample(1:81, 3, replace=F)
#gives: 12, 39, 27
h$SPECIES[12]
h$SPECIES[39]
h$SPECIES[27]


### gotta to be a simplier way to do this... 
#This creates 0 counts for 3 unique data sets based on each species number


COUNT<-rep(0,24)
hsites<-unique(h[c("LON_DEG", "LAT_DEG")])

test<-cbind(hsites$LON_DEG, hsites$LAT_DEG, sites, COUNT)
test<-as.data.frame(test)
names(test)<-c("LON_DEG", "LAT_DEG","EnvSiteID", "COUNT")
test$EnvSiteID<-as.integer(test$EnvSiteID)


        ### todo:  Turn below into loop ###

#Aceteon melampoides
h.12<-h[which(h$SPECIES== h$SPECIES[12]) , ]
h.12<-h.12[, c(24,23,1,20)]
b<-setdiff(test[,(1:3)], h.12[,(1:3)])
b$COUNT<-(rep(0,length(b$LON_DEG)))
h.12<-union(b, h.12)

#Benthonella tenella
h.39<-h[which(h$SPECIES== h$SPECIES[39]) , ]
h.39<-h.39[, c(24,23,1,20)]
b<-setdiff(test[,(1:3)], h.39[,(1:3)])
b$COUNT<-(rep(0,length(b$LON_DEG)))
h.39<-union(b, h.39)

#Theta lyronuclea
h.27<-h[which(h$SPECIES== h$SPECIES[27]) , ]
h.27<-h.27[, c(24,23,1,20)]
b<-setdiff(test[,(1:3)], h.27[,(1:3)])
b$COUNT<-(rep(0,length(b$LON_DEG)))
h.27<-union(b, h.27)

#############################   Clear Varaibles:    #############################

rm(hwd, hpts,b,hloc,h.utm, hsites, nw, o, oloc, osites, test,h, COUNT, res, sites, opts,d, NE, neopts, NEbound)

############################  Non-Spatial Analysis    ###################################
b<-h.div$Date
b<-na.omit(b)

div_by_dept<-lm(h.div$Diversity ~ h.div$Depth)

plot(h.div$Diversity ~ h.div$Depthm, xlab=" Depth", ylab="Number of Species", main="Diversity")
abline(div_by_dept)

############################  Geostats NE   ###################################

############################    Diversity   ###################################
d<-as.geodata(h.div, coords.col =3:4, data.col = 2 )
d.var<-(d$data)
d.v <- variog(d)
d.v4 <- variog4(d)

wls.a<-variofit(a.v, ini= c(50,4*10^5), cov.model="matern", nugget=25)
wlspe.a<- variofit(a.v, ini= c(250,4*10^5), cov.model="powered.exponential", nugget=25)
reml.a <- likfit(a.v, ini = c(30,2*10^5), nugget = 25, cov.model="exponential",lik.method = "RML")

#############################   Abundances   ###################################



          #species a
a<-as.geodata(h.12, coords.col=1:2,data.col=4)
a.var<-(a$data)
a.v <- variog(a)
a.v4 <- variog4(a)

wlsm.a <- variofit(a.v, ini= c(35,2*10^5), cov.model="matern", nugget=20)
reml.a <- likfit(a, ini.cov.pars = c(30,2*10^5), nug = 20, cov.model="exponential",lik.method = "RML")

summary(wlsm.a)


#species b
b<-as.geodata(h.39, coords.col=1:2,data.col=4)
b.var<-(b$data)
b.v<-variog(b)
b.v4 <- variog4(b)
wlsm.b <- variofit(b.v, ini= c(600,4*10^5), cov.model="matern", nugget=500)

reml.b <- likfit(b, ini.cov.pars = c(30,2*10^5), nug = 20, cov.model="exponential",lik.method = "RML")

summary(wlsm.b)


3.693605e+05/1000

#species T

t<-as.geodata(h.27, coords.col=1:2,data.col=4)
summary(t)
t.var<-(t$data)
t.v<-variog(t)
t.v4 <- variog4(t)
wlsm.t <- variofit(t.v, ini= c(0.3,5*10^5), cov.model="matern", nugget=1.7)


summary(wlsm.t)

#############################   Krigging   ###################################
                  ####     3 Species Abundance        ###


##set area to predict... to big.. need to adjust
xr<-max(h.12$LON_DEG)-min(h.12$LON_DEG)
yr<-max(h.12$LAT_DEG)-min(h.12$LAT_DEG)

xo<-(min(h.12$LON_DEG)-2000)
xn<-(max(h.12$LON_DEG)+2000)
yo<-(min(h.12$LAT_DEG)-2000)
yn<-(max(h.12$LAT_DEG)+2000)

x<-seq(xo,xn,100)
y<-seq(yo,yn,100)

d1 <- expand.grid(x=x,y=y)
a.conv <- krige.conv(a, loc=d1, krige=krige.control(obj.m = wlsm.a))
plot(a.conv,loc=d1,values=a.conv$predict,
      main="Conventional kriging")

#### stopped here ###
                                ####
                                ####

#############################   PLOTS   ###################################

setwd(plotsave)
dev.off()
par(mfrow=c(1,1))

#############################   Tables   ###################################

# table 1

#write.table(h.div, "table1.csv", row.names = F)

#############################   2d Mapping   ###################################
#color ramp for 2d image
colorz<-c("#2D3184", "#254289", "#1C518F", "#115F96", "#046C9C", "#0078A2", "#0084A7", "#068FAB", "#189AAF", "#28A4B3", "#45B6B8", "#54BEBA", "#62C5BC", "#70CCBD", "#7DD2BF", "#8BD8C0", "#97DDC1", "#A3E2C3", "#AFE5C4", "#BAE9C6", "#C4EBC8", "#CDEECA", "#D6F0CD", "#DEF1D0", "#E4F2D3", "#EAF2D6", "#EEF2DA", "#F2F2DE", "#F3F1E4")
ramp<-colorz

tiff("NE_30__final_isomap.tiff", res=600, height=4800, width=4800)
plot(NE, bpal=ramp, main =" Figure 1: NW Atlantic (New England) 24 sites. 200m isobaths", image=T, deep=c(-6500,0), shallow=c(-50,0), step=c(200,0), 
     lwd=c(0.4,0.8), lty=c(1,1))
scaleBathy(NE, deg=2, y=42.5,x=-75)
points(neopts, pch=21, col="orange",bg=col2alpha("red",.4),cex=.8)
dev.off()

#### sites between 300- 4000
shelf<-d[d$Depth> 3000 & d$Depth< 4000,]
ssites<-cbind(shelf$Longitude,shelf$Latitude)
spts<-SpatialPoints(ssites)
proj4string(spts)<-CRS("+proj=longlat +datum=WGS84")



tiff("NE_30__isomap_break.tiff", res=300, height=2400, width=2400)
plot(NE, bpal=ramp, main =" Figure 6: 6 sites at Continetal Rise. 200m isobaths", image=T, deep=c(-6500,0), shallow=c(-50,0), step=c(200,0), 
     lwd=c(0.4,0.8), lty=c(1,1))
scaleBathy(NE, deg=2, y=42.5,x=-75)
points(spts, pch=21, col="orange",bg=col2alpha("red",.4),cex=1.4)
dev.off()

#############################   3d lattice   ###################################
              ### hack to label axis correctly..... fix this


#### wireframe color bar: ##

shade.col.fun <- trellis.par.get("shade.colors")$palette
shade.colors <- shade.col.fun(0.5, 0.5, seq(0, 1, length = 100))



Depth<-unclass(NE)
tiff("NE30_TOP.tiff", res=300, height=2400, width=2400)
wireframe(Depth, scales=list(y=list(rot=90)),  xlab = "Longitude", ylab = "Latitude", shade=T, aspect = c(1,1), screen = FALSE,
          colorkey = list(col = shade.colors,at = do.breaks(range(Depth), 100)),
          par.settings = list(axis.line = list(col = "transparent")),
          par.box = c(col = rgb(0,0,0,0,0.2)))
dev.off()
rm(Depth)

Latitude<-unclass(NE)
tiff("NE30_side.tiff", res=300, height=2400, width=2400)
wireframe(Latitude, scales=list(y=list(rot=90)), xlab = "Longitude", ylab = "Depth", shade=T, aspect = c(1,1), screen = list(z = 1, x = -55),
          colorkey = list(col = shade.colors,at = do.breaks(range(Latitude), 100)),
          par.settings = list(axis.line = list(col = "transparent")),
          par.box = c(col = rgb(0,0,0,0,0.2)))
dev.off()

rm(Latitude)

#############################   Non-Spatial. Analysis -plots    ###################################



tiff("Site Gear Type.tiff", res=300, height=2400, width=2400)
plot(h$GEAR_NAME, main= "Site Gear Type N= 1")

dev.off()


tiff("Diversity_depth.tiff",  res=300, height=2400, width=2400)
plot(h.div$Diversity ~ h.div$Depth, xlab=" Depth", ylab="Number of Species", main="Diversity",
     pch=25, col='blue', bg='black',cex=1.2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
abline(div_by_dept, col="red")
dev.off()

h.lm<-lm(h$COUNT~h$MAX_DEPTH)
tiff("Counts by Depth", res=300, height=4800, width=4800)
plot(h$COUNT~ sort(h$MAX_DEPTH), xlab= "Depth", ylab="Count", pch=25, bg="blue", cex=2, cex.lab=2, cex.axis=2)
abline(h.lm)
dev.off()

h$MAX_DEPTH<-h$MAX_DEPTH* -1
tiff("depth_range.tiff", res=300, height=2400, width=2400)
plot((unique(h$MAX_DEPTH)),main="Depth Range of Sites", ylab="Depth: meters below sea level", xlab="Sites", pch=25, col='blue', bg='black',cex=2, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()
h$MAX_DEPTH<-h$MAX_DEPTH *-1


tiff("Rarity in Abundance", res=300, height=2400, width=2400)
qplot(h$COUNT, data=h, geom="histogram", binwidth= 1, xlab = "indviduals of one species per site", cex=2, cex.lab=.5, cex.axis=2)
dev.off()

#############################   Krigging-plots- Diversity Variograms   ###################################

tiff("explortatroy_diversity.tiff", res=300, height=2400, width=2400)
plot(d)
dev.off()


tiff("Directional_&_Omni_Diversity.tiff", res=300, height=2400, width=2400)
par(mfrow=c(1,2))
plot(d.v)
abline(h=var(d$data), col="orange")
plot(d.v4)
abline(h=var(d$data), col="orange")
dev.off()

tiff("Omnidirectional_div_mattern.tiff",  res=300, height=2400, width=4800)
par(mfrow=c(1,1))
plot(d.v, pch=19, col="blue", cex=1, cex.lab=1.5, cex.axis=1,  cex.sub=1.5)
abline(h=var(t$data), col="red")
lines(wlsm.a, lwd=1)
dev.off()
par(mfrow=c(1,1))

#############################   Krigging-plots- 3 Species Abundance Variograms   ###################################
                  
tiff("Directional.tiff",  res=300, height=2400, width=4800)
par(mfrow=c(1,3))
plot(t.v4)
abline(h=var(t$data), col="red")  

plot(a.v4)
abline(h=var(a$data), col="orange")  

plot(b.v4)
abline(h=var(b$data), col="blue") 
dev.off()
par(mfrow=c(1,1))

dev.off()

tiff("Omnidirectional.tiff",  res=300, height=2400, width=4800)
par(mfrow=c(1,3))
plot(t.v, pch=19, col="red", cex=2, cex.lab=1.5, cex.axis=1,  cex.sub=1.5)
abline(h=var(t$data), col="red")
lines(wlsm.t, lwd=1)

plot(a.v, pch=24, col="black", bg="orange", cex=2, cex.lab=1.5, cex.axis=1,  cex.sub=1.5)
abline(h=var(a$data), col="orange")
lines(wlsm.a, lwd=1)

plot(b.v, pch=15, col="blue", cex=2, cex.lab=1.5, cex.axis=1,  cex.sub=1.5)
abline(h=var(b$data), col="blue")
lines(wlsm.b, lwd=1)
dev.off()
par(mfrow=c(1,1))



tiff("am_ex.tiff",  res=300, height=2400, width=2400)
plot(a)
dev.off()

tiff("bt_ex.tiff", res=300, height=2400, width=2400)
plot(b)
dev.off()

tiff("Tl_ex.tiff", res=300, height=2400, width=2400)
plot(t)
dev.off()

############################    REFERENCES ####################################




# color ramp-wire:  https://stat.ethz.ch/pipermail/r-help/2007-October/143136.html


# points: https://stat.ethz.ch/pipermail/r-help/2005-November/083086.html
# http://r.789695.n4.nabble.com/lattice-wireframe-quot-eats-up-quot-points-how-to-make-points-on-wireframe-visible-td3334836.html

# show transect of sidescan:

### https://cran.r-project.org/web/packages/marmap/vignettes/marmap-DataAnalysis.pdf

# http://robinlovelace.net/r/2014/07/29/clipping-with-r.html

# http://www2.hawaii.edu/~bmgenco/r.html

# http://spatialreference.org/ref/epsg/32619/


