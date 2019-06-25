###Heatmap Generation using Generic Mapping Tools (GMT)
#v0.2
#by: Livia Tran 

#required libraries 
require(stringr)
require(data.table)
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





############HEATMAP GENERATION SCRIPT  
require(gmt)
require(DescTools)

#puts motif_coords long, lat, and allele frequency columns into a new variable to be converted
#into gmt format (xyz file)
gmt_converted_data<-motif_coords[,c(1,2,3,4)]

#rerranges format 
gmt_converted_data<-gmt_converted_data[,c("long", "lat", "allele_freq", "pop")]

#calls function to convert R object gmt_converted_data to GMT formatted data 
r2gmt(gmt_converted_data, "DRB1*26F~28E~30Y.xyz")

#finds blockmean of data 
#-60 for lower bound lat
#+80 for upper bound lat
#I for every 3x3 mean value, as obtained from gmt.sh file
gmt.system("blockmean DRB1*26F~28E~30Y.xyz -R-180/180/-60/80 -I3 > DRB1*26F~28E~30Y.block")

#grids table using surface 
gmt.system("surface DRB1*26F~28E~30Y.block -R-180/180/-60/80 -Gmotif.grd -I1.8 -T.7 -Ll0")

#creates basemap with correct motif 
gmt.system("psbasemap -JM6i -R-180/180/-60/80 -B0:.`echo DRB1*26F~28E~30Y.xyz| sed 's/.xyz//g'`: -K > basemap.ps")

#overlays relevant continents onto basemap 
gmt.system("pscoast -JM6i -R-180/180/-60/80 -A30000 -B0 -G200 -W0.25p -O -K >> basemap.ps")	

#gets upperbound for allele frequencies
gmt.system("awk '{print $3}' DRB1*26F~28E~30Y.xyz | sort -r | head -1 > upperbound")

#makes a vector called cpt_interval which contains max allelic information and decile increment needed
#to form deciles
#uses readLines to obtain upperbound information from bash, and rounds it to the nearest 0.5 
#uses readLines to obtain decile needed
cpt_interval<-c(RoundTo(as.numeric(readLines("upperbound")), 0.5), RoundTo(as.numeric(readLines("upperbound")), 0.5)/10)

#creates a vector called decile_interval, which gives decile increments based on cpt_interval information
decile_interval<-seq(0, cpt_interval[1], cpt_interval[2])

#writes max allelic frequency (rounded) and increment needed to form deciles 
write(cpt_interval, "max_cpt")

#writes decile interval without any line breaks
cat(decile_interval, file="deciles")

#makes custom CPT with increments of max_frequency/10 for deciles
#specifically calls on max_cpt for max frequency and decile increments 
gmt.system("makecpt -Cseis -Iz -T0/`awk '{print $1}' max_cpt`/`awk '{print $2}' max_cpt` > decile.cpt")

#adds color scale to basemap based on cpt provided
gmt.system("psscale -D0.1i/1.1i/2i/0.3i -Cdecile.cpt -Np -L -O -K >> basemap.ps")

#overlays more coastlines with pscoast
gmt.system("pscoast -JM6i -R-180/180/-60/80 -A100000 -Gc -O -K >> basemap.ps")

#clips/masks map areas with no data table coverage -- radius of influence increased to 900 km 
gmt.system("psmask DRB1*26F~28E~30Y.xyz -R-180/180/-60/80 -I3 -JM6i -S900k -O -K >> basemap.ps")

#grids .grd file onto map 
gmt.system("grdimage motif.grd -Cdecile.cpt  -JM6i -R-180/180/-60/80 -O -K >> basemap.ps")

#makes a contour map using a .grd file 
gmt.system("grdcontour motif.grd -JM6i -Cdecile.cpt -A- -R-180/180/-60/80 -O -K >> basemap.ps") 

#plots longtitude/latitude coordinates onto basemap
gmt.system("psxy DRB1*26F~28E~30Y.xyz -R-180/180/-60/80 -JM6i -A -G255 -W0.5p -Sc.05 -O -K >> basemap.ps")

#overlays coastlines with pcoast 
gmt.system("pscoast -JM6i -R-180/180/-60/80 -A10000 -W.5p -O  >> basemap.ps")



