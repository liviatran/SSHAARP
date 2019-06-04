###heatmap generation script
#v0.1
#by: Livia Tran 

###conversions for latitude and longitude
#converts coordinates to appropriate enumerations based on cardinal directions -- inserts back
#into appropriate rows in HMD
#NORTH and EAST are positive
  #N = above equator
  #E = east the prime meridian 
#SOUTH AND WEST ARE NEGATIVE
  #S = below equator
  #W = west of the prime meridian 

#latitude conversions
if(any((grepl("S", HMD$latit))==TRUE)){
HMD$latit[which((grepl("S", HMD$latit))==TRUE)]<-as.numeric(paste("-", str_extract(HMD$latit[which((grepl("S", HMD$latit))==TRUE)],"\\-*\\d+\\.*\\d*"), sep=""))}

if(any((grepl("N", HMD$latit))==TRUE)){
HMD$latit[which((grepl("N", HMD$latit))==TRUE)]<-as.numeric(str_extract(HMD$latit[which((grepl("N", HMD$latit))==TRUE)],"\\-*\\d+\\.*\\d*"))
}

#longitude conversions
if(any((grepl("W", HMD$longit))==TRUE)){
  HMD$longit[which((grepl("W", HMD$longit))==TRUE)]<-as.numeric(paste("-", str_extract(HMD$longit[which((grepl("W", HMD$longit))==TRUE)],"\\-*\\d+\\.*\\d*"), sep=""))}

if(any((grepl("E", HMD$longit))==TRUE)){
  HMD$longit[which((grepl("E", HMD$longit))==TRUE)]<-as.numeric(str_extract(HMD$longit[which((grepl("E", HMD$longit))==TRUE)],"\\-*\\d+\\.*\\d*"))
}

#filters out admixed, migrant OTH populations
HMD<-HMD[HMD$contin!="OTH",]

coords<-sapply(HMD$popname, function(x) list())

for(i in 1:length(coords)){
coords[[i]]<-HMD[,c(3,4,5)][HMD[,1] %in% names(coords)[[i]],]}


motif_coords<-melt(coords, id.vars=c("latit", "longit"))


motif_coords$value<-as.numeric(as.character(motif_coords$value))
motif_coords$latit<-round(as.numeric(motif_coords$latit), digits = 1)
motif_coords$longit<-round(as.numeric(motif_coords$longit), digits=1)


colnames(world)<-c("longit", "latit", "group", "order", "region", "subregion")

############heatmap generation 
library(data.table)
library(dplyr)
library(ggplot2)
library(grid)
library(reshape)
library(wesanderson)

##gets world map data
world <- map_data('world') %>% data.table()
world <- world[region!='Antarctica',]


#plots blank world map 
 
  g <- ggplot(world, aes(long, lat)) + 
  geom_polygon(aes(group=group),fill="white",colour="black",size=0.1) +
  coord_equal() + 
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  labs(x='Longitude', y='Latitude') +
  theme_bw()


g

world$long<-round(world$long, digits=1)
world$lat<-round(world$lat, digits=1)

  

  
  
  
  

#map with color bar as legend
g+
  geom_point(data=motif_coords, aes(x = longit, y = latit, color = factor(motif_coords$value)), size=1.5)



  scale_colour_gradientn(colours = wes_palette("FantasticFox1"),
                         guide="colorbar", limits = c(0, round(max(motif_coords$value), digits=3)))
                         
                         
                         
                         
                         
                         #labels = c(0,round(max(motif_coords$value)/10:max(motif_coords$value), digits=3)),
                        #breaks= c(0,round(max(motif_coords$value)/10:max(motif_coords$value), digits=3)))+
          theme(legend.key.size = unit(3.0, "cm"))
      

          
         
g+
  geom_tile(data=motif_coords, aes(x=longit, y=latit, color=motif_coords$value))

#ggplot(data = motif_coords, aes(x=longit, y=latit, fill=value))+


g + 
  geom_tile(data=motif_coords,aes(longit, latit, fill=value)) +
  scale_fill_gradient(low="blue", high="red", name='Yield Gap (ton/ha)') +
  theme(
    legend.position = 'bottom',
    legend.key.size = unit(1, "cm")
  )

                                                             

    


