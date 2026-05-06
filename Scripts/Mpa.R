# Fig 1 Maps 
library(PBSmapping)
library(ggplot2)
library(maps)
library(mapdata)
library(mapproj)
library(ggspatial)



GW_Map<- read.csv(here("./data/GW_Biol.csv")) 

colnames(GW_Map) <- c("vial", "MMSNW" , "ID","BlubberDEPTH","Age.class","UME","Yr","Mo","Da","Month",
                   "Latitude.N","Longitude.W","State", "CarCode" ,  "BodyCond","Adjusted_BC",  "Sex", "COD","Lipid",      
                   "Cts", "Cts CV", "Ctn" ,  "Ctn CV",  "Aldo", "Aldo CV")


GW_Map <- GW_Map %>%
  select(MMSNW, ID, Latitude.N, Longitude.W, State, Sex)
GW_Map$Sex = as.factor(GW_Map$Sex)
(rangeLAT<- range(GW_Map$Latitude.N, na.rm=TRUE)) #46.33 60.82
(rangeLong<- range(GW_Map$Longitude.W, na.rm=TRUE)) #-161.76 -122.26
summary(GW_Map$Sex) #22 females and 13 males 

GW_Map$State = as.factor(GW_Map$State)

#fig 1: PNW
xlim = c(-170,-110)
xmin=-170
xmax=-110
ymin=45
ymax=65
ylim = c(45, 65)
par(mar=c(4, 7.0, 4, 8.0))  
worldmap = map_data("world")
#  Change the names of the worldmap data frame to be consistent with the expected input for the clipPolys function.
names(worldmap) <- c("X","Y","PID","POS","region","subregion")

# just the range you want to plot. Otherwise, it has a tough time connecting the polygon lines.
worldmap = clipPolys(worldmap, xlim=xlim,ylim=ylim, keepExtra=TRUE)
statemap = map_data("state")
names(statemap) <- c("X","Y","PID","POS","region","subregion")
statemap = clipPolys(statemap, xlim=xlim,ylim=ylim, keepExtra=TRUE)  

??scalebar
Map <-ggplot() +
  coord_map(xlim=xlim,ylim=ylim) +
  geom_polygon(data=worldmap,aes(X,Y,group=PID),
               fill = "gray98",color="black", lwd=0.5)+
  geom_polygon(data=statemap,aes(X,Y,group=PID),
               fill = "gray98",color="black", linetype="twodash")+
  xlab(expression(paste(Longitude^o,~'W'))) +
  ylab(expression(paste(Latitude^o,~'N'))) +
  # scalebar(location="topright", transform=T, dist_unit = "km", dist=500,
  #          x.min=xmin+30, x.max=xmax, y.min=ymin+3, y.max=ymax, st.color="black",
  #          st.size=3) +
  annotation_north_arrow(location ="bl")+
  theme_bw() +theme(legend.position = "none")+
  theme(axis.title.y=element_text(size=18)) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.text.x=element_text(size=16)) +
  theme(axis.text.y=element_text(size=16)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) 

Map

Fig1=Map +
  annotate("text", x = -150, y =63, label = "AK", col="black", size=4, fontface= 2)+
  annotate("text", x = -120, y =47, label = "WA", col="black", size=4, fontface= 2)

Fig1
ggsave(here("./output/Map.png"),Fig1, width=6, height=6, dpi=300)


summary(GW_Map$State)
# 7 alaska and 28 washington

# figure 1A Washington state
xlim = c(-128,-120)
xmin=-128
xmax=-120
ymin=46
ymax=50
ylim = c(46, 50)
par(mar=c(4, 7.0, 4, 8.0))  
worldmap = map_data("world")
#  Change the names of the worldmap data frame to be consistent with the expected input for the clipPolys function.
names(worldmap) <- c("X","Y","PID","POS","region","subregion")

# just the range you want to plot. Otherwise, it has a tough time connecting the polygon lines.
worldmap = clipPolys(worldmap, xlim=xlim,ylim=ylim, keepExtra=TRUE)
statemap = map_data("state")
names(statemap) <- c("X","Y","PID","POS","region","subregion")
statemap = clipPolys(statemap, xlim=xlim,ylim=ylim, keepExtra=TRUE)  

Fig1a_WA <- ggplot() +
  coord_map(xlim=xlim,ylim=ylim) +
  geom_polygon(data=worldmap,aes(X,Y,group=PID),
               fill = "gray98",color="gray29", lwd=0.5)+
  xlab(expression(paste(Longitude^o,~'W'))) +
  ylab(expression(paste(Latitude^o,~'N'))) +
  # scalebar(location="topright", transform=T, dist_unit = "km", dist=500,
  #          x.min=xmin+30, x.max=xmax, y.min=ymin+3, y.max=ymax, st.color="black",
  #          st.size=3) +
  annotation_north_arrow(location ="bl")+
  geom_point(data=subset(GW_Map, GW_Map$State=="WA"), aes(Longitude.W, Latitude.N, shape=Sex, fill=Sex), size=2)+
  scale_shape_manual(values=c(21,25))+
  scale_fill_manual(values= c("black", "red"))+
  theme_bw() +theme(legend.position = "none")+
  theme(axis.title.y=element_text(size=18)) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.text.x=element_text(size=16)) +
  theme(axis.text.y=element_text(size=16)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) 
Fig1a_WA

# ggsave(here("./output/MapWA.png"),Fig1a_WA, width=6, height=6, dpi=300)

#Alaska map

xlim = c(-170,-130)
xmin=-130
xmax=-170
ymin=50
ymax=65
ylim = c(50,65)
par(mar=c(4, 7.0, 4, 8.0))  
worldmap = map_data("world")
#  Change the names of the worldmap data frame to be consistent with the expected input for the clipPolys function.
names(worldmap) <- c("X","Y","PID","POS","region","subregion")

# just the range you want to plot. Otherwise, it has a tough time connecting the polygon lines.
worldmap = clipPolys(worldmap, xlim=xlim,ylim=ylim, keepExtra=TRUE)
statemap = map_data("state")
names(statemap) <- c("X","Y","PID","POS","region","subregion")
statemap = clipPolys(statemap, xlim=xlim,ylim=ylim, keepExtra=TRUE)  

Fig1a_AK <- ggplot() +
  coord_map(xlim=xlim,ylim=ylim) +
  geom_polygon(data=worldmap,aes(X,Y,group=PID),
               fill = "gray98",color="gray29", lwd=0.5)+
  xlab(expression(paste(Longitude^o,~'W'))) +
  ylab(expression(paste(Latitude^o,~'N'))) +
  # scalebar(location="topright", transform=T, dist_unit = "km", dist=500,
  #          x.min=xmin+30, x.max=xmax, y.min=ymin+3, y.max=ymax, st.color="black",
  #          st.size=3) +
  annotation_north_arrow(location ="bl")+
  geom_point(data=subset(GW_Map, GW_Map$State=="AK"), aes(Longitude.W, Latitude.N, shape=Sex, fill=Sex), size=2)+
  scale_shape_manual(values=c(21,25))+
  scale_fill_manual(values= c("black", "red"))+
  theme_bw() +theme(legend.position = "none")+
  theme(axis.title.y=element_text(size=18)) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.text.x=element_text(size=16)) +
  theme(axis.text.y=element_text(size=16)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) 
Fig1a_AK
# ggsave(here("./output/MapAK.png"),Fig1a_AK, width=6, height=6, dpi=300)
Fig1_map<- ggarrange( Fig1a_AK, Fig1a_WA,nrow= 1, ncol = 2)

ggsave(here("./output/Fig1_map.png"),Fig1_map, width=8, height=6, dpi=600)
