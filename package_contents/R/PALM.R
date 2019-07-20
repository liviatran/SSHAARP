#'Population Allele Locating Mapmaker
#'
#'Produces a heatmap from the manipulated Solberg dataset by using the Generic Mapping Tools (GMT) R Package, which is an interface between R and the GMT Map-Making software. GMT commands are called via a Bash script through using the gmt.system() function in the GMT R package.
#'
#'@param gdataset Solberg dataset. Can be found at http://pypop.org/popdata/ under results.zip.
#'@param motif An amino acid motif in the following format: Locus*AA_1*AA_2*AA_3.
#'@param color Should the heatmap produced be in color? The default for color is set to true for a color heatmap -- for a map in grayscale, set color==FALSE.
#'@param filter_migrant Should OTH and migrant populations be filtered out of the dataset? The default for filter_migrant is set to true -- set this parameter equal to FALSE to keep migrant and OTH populations in the dataset.
#'
#'@importFrom gmt gmt.system
#'@importFrom DescTools RoundTo
#'
#'@return A .jpg file of the produced heatmap is output to the user's working directory.
#'
PALM<-function(gdataset, motif, color=TRUE, filter_migrant=TRUE){

  #uses dataSubset to read and manipulate the Solberg dataset
  solberg_DS<-dataSubset(gdataset, motif)

  #makes an empty list named unique_AWM, where the name of each element is after a unique AWM,
  #which is acquired by using the motif_finder function
  unique_AWM<-sapply(unique(motif_finder(motif)$trimmed_allele), function(x) NULL)

  #finds unique_AWMs in Solberg dataset and extracts the allele frequency and locus_allele column
  for(y in 1:length(unique_AWM)){
    unique_AWM[[y]]<-solberg_DS[,c(2,3,4,5,6,11,14)][solberg_DS$locus_allele %in% names(unique_AWM[y]),]}

  #subsets out locus_allele pairs with the motif from the alignment but aren't present in the Solberg ds
  unique_AWM<-unique_AWM[sapply(unique_AWM, nrow)>0]

  #melts dataframes in list into one big dataframe
  unique_AWM<-melt(unique_AWM, id.vars=c("popname", "contin", "complex", "latit", "longit", "allele.freq", "locus_allele"))

  #reorders dataframe based on popname alphabetical order
  unique_AWM<-unique_AWM[order(unique_AWM$popname),]

  #makes row names NULL so they are back in sequential enumeration
  row.names(unique_AWM)<-NULL

  #creates a variable named Population Allele Frequencies (PAF), where each element is named after a
  #unique popname for the targeted locus
  #makes each element a list to take in allele frequencies in the next for loop
  PAF<-sapply(unique(solberg_DS$popname), function(x) NULL)

  #for loop for finding allele frequencies for each named element in PAF
  #if multiple entries are present for a given population name, the allele frequencies are added up
  #if the population does not have the motif, the allele frequency is 0, since the entry is not found in unique_AWM
  for(i in 1:length(PAF)){
    PAF[[i]]<-sum(as.numeric(unique_AWM$allele.freq[which((grepl(names(PAF[i]), unique_AWM$popname))==TRUE)]))}

  #melts PAF into a two columned df
  PAF<-melt(PAF)

  #renames column names
  colnames(PAF)<-c("allele_freq", "popname")

  #merges PAF with solberg_DS data to a new variable, tbm_ds
  #which contains summed up allele frequencies, and relevant contin, complex, locus*allele, and coordinate information
  tbm_ds<-merge(PAF, solberg_DS[!duplicated(solberg_DS$popname),], by="popname")[,c("popname","contin", "complex", "latit", "longit", "allele_freq")]

  #converts coordinates to proper enumerations
  tbm_ds<-coordinate_converter(tbm_ds)

  #filters out migrant populations if filter_migrant==TRUE
  if(filter_migrant==TRUE){
    #filters out admixed, migrant OTH populations
    motif_map_df<-tbm_ds[tbm_ds$contin!="OTH",]

    #filters out migrant populations in complexity column
    motif_map_df<-motif_map_df[which((grepl("mig", motif_map_df$complex))==FALSE),]}

  #specifies certain rows from motif_map_df to go into a new variable, gmt_converted_data
  gmt_converted_data<-motif_map_df[,c(1,4,5,6)]

  #rerranges format
  gmt_converted_data<-gmt_converted_data[,c("longit", "latit", "allele_freq", "popname")]

  #calls function to convert R object gmt_converted_data to GMT formatted data
  #outputs converted data into local environment
  r2gmt(gmt_converted_data, "motif.xyz")

  #writes motif to working directory so bash script can extract
  #and use it as the title for the map
  write(motif, "motif")

  #finds blockmean of data
  #-60 for lower bound lat
  #+80 for upper bound lat
  #I for every 3x3 mean value, as obtained from gmt.sh file
  gmt.system("blockmean motif.xyz -R-180/180/-60/80 -I3 > motif.block")

  #grids table using surface
  gmt.system("surface motif.block -R-180/180/-60/80 -Gmotif.grd -I0.5 -T.7 -Ll0")

  #creates basemap with correct motif
  gmt.system("psbasemap -JM6i -R-180/180/-60/80 -B0:.`cat motif`: -K > basemap.ps")

  #uses white background to fill landmasses if color=T
  if(color==TRUE){
    #overlays relevant continents onto basemap
    gmt.system("pscoast -JM6i -R-180/180/-60/80 -A30000 -B0 -G200 -W0.25p -O -K >> basemap.ps")}

  #uses a hashed background to fill landmasses if color=F
  if(color==FALSE){
    gmt.system("pscoast -JM6i -R-180/180/-60/80 -A30000 -B0 -Gp61 -W0.25p -O -K >> basemap.ps")
  }

  #gets upperbound for allele frequencies
  gmt.system("awk '{print $3}' motif.xyz | sort -r | head -1 > upperbound")

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

  #uses seis color palette if color=T
  if(color==TRUE){
    #makes custom CPT with increments of max_frequency/10 for deciles
    #specifically calls on max_cpt for max frequency and decile increments
    gmt.system("makecpt -Cseis -Iz -T0/`awk '{print $1}' max_cpt`/`awk '{print $2}' max_cpt` > decile.cpt")
  }

  #uses grayscale color palette if color=F
  if(color==FALSE){
    gmt.system("makecpt -Cgray -Iz -T0/`awk '{print $1}' max_cpt`/`awk '{print $2}' max_cpt` > decile.cpt")

  }
  #adds color scale to basemap based on cpt provided
  gmt.system("psscale -D0.1i/1.1i/2i/0.3i -Cdecile.cpt -Np -L -O -K >> basemap.ps")

  #overlays more coastlines with pscoast
  gmt.system("pscoast -JM6i -R-180/180/-60/80 -A100000 -Gc -O -K >> basemap.ps")

  #clips/masks map areas with no data table coverage -- radius of influence increased to 900 km
  gmt.system("psmask motif.xyz -R-180/180/-60/80 -I3 -JM6i -S850k -O -K >> basemap.ps")

  #grids .grd file onto map
  gmt.system("grdimage motif.grd -Cdecile.cpt  -JM6i -R-180/180/-60/80 -O -K >> basemap.ps")

  #makes a contour map using a .grd file
  gmt.system("grdcontour motif.grd -JM6i -Cdecile.cpt -A- -R-180/180/-60/80 -O -K >> basemap.ps")

  #plots longtitude/latitude coordinates onto basemap
  gmt.system("psxy motif.xyz -R-180/180/-60/80 -JM6i -A -G255 -W0.5p -Sc.05 -O -K >> basemap.ps")

  #calls psmask again to terminate clip path with -C parameter
  gmt.system("psmask motif.xyz -R-180/180/-60/80 -I3 -JM6i -S850k -C -O -K >> basemap.ps")

  #calls pcoast again to re-establish coastlines and -Q parameter to quit clipping
  gmt.system("pscoast -JM6i -R-180/180/-60/80 -A10000 -W0.5 -O -Q  >>  basemap.ps")

  #converts ps map to jpg -- saves into local environment
  #requires Ghostscript in order to execute command
  gmt.system("psconvert basemap.ps -A -Tj -P -Qg4 -E2000")
}
