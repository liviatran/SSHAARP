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

#filters out migrant populations in complexity column 
HMD<-HMD[which((grepl("mig", HMD$complex))==FALSE),]

coords<-sapply(HMD$popname, function(x) list())

#gets coordinate and allele frequency information for each allele
for(i in 1:length(coords)){
coords[[i]]<-HMD[,c(4,5,6)][HMD[,1] %in% names(coords)[[i]],]}

#melts all data together into one dataframe 
motif_coords<-melt(coords, id.vars=c("latit", "longit"))

#gets rid of melt joint variable
motif_coords$variable<-NULL

#sets allele freq to numeric class
motif_coords$value<-as.numeric(as.character(motif_coords$value))

#sets coordinates to numeric class and rounds to 4 digits 
motif_coords$latit<-round(as.numeric(motif_coords$latit), digits = 4)
motif_coords$longit<-round(as.numeric(motif_coords$longit), digits=4)

#renames column names
colnames(motif_coords)<-c("lat", "long", "allele_freq", "pop")


############heatmap generation 
library(gmt)


#adds a size column to determine how big to make plot points
motif_coords$size<-rep("3", nrow(motif_coords))
motif_coords$color<-NULL
#rearranges columns to fit gmt format 
motif_coords<-motif_coords[,c("long", "lat", "size", "pop", "allele_freq")]

#world map generation
pscoast("-R-180/180/-80/80 -JM6i -B30f5/20f2 -Ggrey -K", "heatmap.ps")

#plots population coordinate points onto map 
psxy(motif_coords, cmd="-J -R -Sdp -W1.5p,black -O -K", file="heatmap.ps")

