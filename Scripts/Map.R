# Fig 1 Maps 
library(sf) 
library(ggplot2) 
library(rnaturalearth) 
library(rnaturalearthdata) 
library(ggspatial) 
library(cowplot) 

library(PBSmapping)
library(maps)
library(mapdata)
library(mapproj)


suppressMessages(library(plyr))
suppressMessages(library(here))
suppressMessages(library(dplyr))

##### creating a big map of the west coast
# Load Canada provinces 
canada <- ne_states(country = "canada", returnclass = "sf") 
# Filter British Columbia 
bc <- canada[canada$name == "British Columbia", ] 
# Load land polygons for context 
land <- ne_countries(scale = "medium", returnclass = "sf")
# Change to WA box and AK bo
### AK
ak_lon_min <- -170
ak_lon_max <- -130
ak_lat_min <- 55
ak_lat_max <- 65
wa_lon_min <- -125
wa_lon_max <- -120
wa_lat_min <- 45
wa_lat_max <- 50

ak_box <- st_as_sfc(st_bbox(c(xmin = ak_lon_min, xmax = ak_lon_max, 
                              ymin = ak_lat_min, ymax = ak_lat_max), 
                            crs = 4326))

wa_box <- st_as_sfc(st_bbox(c(xmin = wa_lon_min, xmax = wa_lon_max, 
                              ymin = wa_lat_min, ymax = wa_lat_max), 
                            crs = 4326))
worldmap = map_data("world")
#  Change the names of the worldmap data frame to be consistent with the expected input for the clipPolys function.
names(worldmap) <- c("X","Y","PID","POS","region","subregion")
world_sf <- st_as_sf(worldmap, coords = c("X", "Y"), crs = 4326) 
world_sf <- world_sf |> group_by(PID) |> summarise(geometry = st_combine(geometry)) |> st_cast("POLYGON") 

# Plot 
West_coast<- ggplot() + geom_sf(data = world_sf, fill = "gray92", color = "black", lwd = 0.5) + 
  #geom_sf(data = bc, fill = "white", color = "black") + 
  geom_sf(data = ak_box, fill = "#FDAE6B", color = "black", linewidth = 0.5, alpha = 0.3) +
  geom_sf(data = wa_box, fill = "#8A8FB4", color = "black", linewidth = 0.5, alpha = 0.5) +
  coord_sf(xlim = c(-170, -115), ylim = c(35, 75)) +
  theme_bw()  + 
  xlab(expression(paste(Longitude^o,~'W'))) + 
  ylab(expression(paste(Latitude^o,~'N')))+
  annotation_north_arrow( location = "tl", which_north = "true", style =     north_arrow_fancy_orienteering )+
  annotation_scale( location = "bl", # bottom-left 
   width_hint = 0.25, # relative width 
   text_cex = 1 )+
  annotate("text", x = -120, y = 55, label = "Canada", size = 5, fontface = "bold") + 
  annotate("text", x = -118, y = 45, label = "USA", size = 5, fontface = "bold") 
West_coast

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

# register_stadiamaps("a3cb2891-9920-489c-bc4b-d6f8c5a77866")
world <- ne_countries(scale = "large", returnclass = "sf")
north_america <- world[world$iso_a2 %in% c("CA", "US"), ]
## Get coordinates for each region 

### AK
ak_lon_min <- -170
ak_lon_max <- -130
ak_lat_min <- 55
ak_lat_max <- 65

ak_box <- st_as_sfc(st_bbox(c(xmin = ak_lon_min, xmax = ak_lon_max, 
                                ymin = ak_lat_min, ymax = ak_lat_max), 
                              crs = 4326))

ak_biop <- GW_Map %>% 
  filter(State == "AK") 
ak_shapes = c("F" = 24, "M" = 21)
ak_cols = c("F" = '#999999', "M" = '#E69F00')
ak_plot <- ggplot() +
  geom_sf(data = north_america, fill = "gray95", color = "black") + 
  coord_sf(xlim = c(ak_lon_min, ak_lon_max), ylim = c(ak_lat_min, ak_lat_max), expand = FALSE) +
  geom_point(data = ak_biop, aes(Longitude.W, Latitude.N,  shape= Sex, fill= Sex), size = 2,
             position = position_jitter(width = 0.2, height = 0.2) ) +
  scale_shape_manual(values = ak_shapes, name= "Sex")+
  scale_fill_manual(values = ak_cols, name= "Sex")+
  xlab(expression(paste(Longitude^o,~'W'))) + 
  ylab(expression(paste(Latitude^o,~'N'))) +
  annotation_scale(location = "br",
                   width_hint = 0.2, # relative width
                   text_cex = 0.8) + 
  theme_bw() +
  theme(legend.position= "inside", legend.position.inside = c(.99,.99), # move legend inside 
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2), 
        legend.title = element_text(size = 10), legend.text = element_text(size = 9, margin = margin(5,5,5,5)))
        
        
ak_plot

wa_biop <- GW_Map %>% 
  filter(State == "WA")

wa_shapes = c("F" = 24, "M" = 21)
wa_cols = c("F" = '#999999', "M" = '#E69F00')
wa_plot <- ggplot() +
  geom_sf(data = north_america, fill = "gray95", color = "black") + 
  coord_sf(xlim = c(wa_lon_min, wa_lon_max), ylim = c(wa_lat_min, wa_lat_max), expand = FALSE) +
  geom_point(data = wa_biop, aes(Longitude.W, Latitude.N,  shape= Sex, fill= Sex), size = 2,
             position = position_jitter(width = 0.2, height = 0.2) ) +
  scale_shape_manual(values = wa_shapes, name= "Sex")+
  scale_fill_manual(values = wa_cols, name= "Sex")+
  xlab(expression(paste(Longitude^o,~'W'))) + 
  ylab(expression(paste(Latitude^o,~'N'))) +
  annotation_scale(location = "br",
                   width_hint = 0.2, # relative width
                   text_cex = 0.8) + 
  theme_bw() +
  theme(legend.position= "inside", legend.position.inside = c(.99,.99), # move legend inside 
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2), 
        legend.title = element_text(size = 10), legend.text = element_text(size = 9, margin = margin(5,5,5,5)))


wa_plot 


ggsave(here("./output/Map/WA_map.png"),wa_plot, width=4, height=8, dpi=600)
ggsave(here("./output/Map/AK_map.png"),ak_plot, width=8, height=4, dpi=600)
ggsave(here("./output/Map/WestCoast.png"),West_coast, width=8, height=8, dpi=600)

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
