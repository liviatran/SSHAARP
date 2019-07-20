#PALM - Popuation Allele Locating Mapmaker script 
#v 1.3
#By: Liv Tran
#7/12/19

#This script contains software developed to get allele frequencies from a 
#population for a given motif from the Solberg dataset, and making heatmaps from 
#allele frequency data by using the GMT R package and calling GMT commands with 
#bash script through gmt.system() in the GMT R package
#output is a jpg heatmap that is saved to the user's working environment 

##REQUIRED PACKAGES 
require(data.table)
require(stringr)
require(gtools)
require(BIGDAWG)
require(gmt)
require(DescTools)

###REQUIRED FUNCTIONS:
#function to count spaces in between regions of interest
#to determine where  start for the alignment sequence 
countSpaces <- function(x){
  counter <- 0
  coll <- numeric()
  vec <- strsplit(x," ")[[1]]
  for(i in 1:length(vec)){
    if (vec[i]==""){
      counter <- counter+1
    }
    else{
      if (counter!=0) coll <- c(coll,counter)
      counter <- 1
    }
  }
  coll
}

######AA_segments_maker
#extracts alignment sequence information for a given locus from the ANHIG/IMGTHLA Github Repo
#to produce a dataframe, where each row is specific to an allele for that locus,
#and each column is the amino acid position, followed by which specific amino acid that allele
#has in that position
#the first 4 columns are locus, allele, trimmed allele, and allele_name 
AA_segments_maker<-function(loci){
#creates empty variables for future for loops
start<-end<-alignment<-list()
  
#creates empty variables where each element is named after the used loci 
  
#empty variables for correspondence table 
pepsplit<-refexon<-AA_aligned<-AA_segments<-inDels<-corr_table<-cols<-downloaded_segments<-w<-alignment_positions<-alignment_length<-alignment_start<-prot_extractions<-refblock_number<-end_char<-space_diff<-sapply(loci, function(x) NULL)
    

#begin for loop   
  for(i in 1:length(loci)){
    #downloads relevant locus alignment file -- readLines allows for space preservation, which is important in
    #finding where the alignment sequence starts 
    alignment[[loci[i]]] <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",paste(ifelse(loci[[i]]=="DRB1","DRB",loci[[i]]),"_prot.txt",sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
    
    #alters alignment file by cutting out non-pertinent information in beginning
    #and endind of alignment file 
    alignment[[loci[i]]] <- head(alignment[[loci[i]]],-3)
    alignment[[loci[i]]] <- tail(alignment[[loci[i]]],-7)
    
    #see countSpaces function at beginning of script 
    #Counts difference between Prot to -30 and beginning of Prot to -30 + 1 due to zero number indexing to find where
    #the alignment sequence actually starts 
    space_diff[[loci[i]]]<-(nchar(strsplit(alignment[[loci[i]]][3], " ")[[1]][2])+countSpaces(alignment[[loci[i]]][3])[2]+1)-countSpaces(alignment[[loci[i]]][2])[1]
    
    #reduces repeated whitespace in alignment file and removes rows with empty values for proper
    #start and stop subsetting 
    alignment[[loci[i]]] <-str_squish(alignment[[loci[i]]])
    alignment[[loci[i]]] <-alignment[[loci[i]]][-which(alignment[[loci[i]]] == "")]
    
    #determines positions of "Prot" and the end of that reference block segment
    start[[loci[i]]]<-as.numeric(grep("Prot", alignment[[loci[i]]]))
    end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,length(alignment[[loci[i]]])))
    
    #counts number of characters in the very last allele to add onto the last Prot enumeration block
    #to obtain end length 
    end_char[[loci[i]]]<-nchar(sapply(strsplit(gsub(" ", "", sub(" ", "~", str_squish(tail(alignment[[loci[i]]], 1)))), "~"), "[", 2))-1
    
    #extracts rows with "Prot" and reference sequence position information 
    #extracts only relevant reference sequence positions
    #NOTE: the first row containing "Prot" contains two numbers -- -30 and 1 -- where only -30, is extracted,
    #as the actual sequence start will always be 1 
    for (j in 1:length(start[[loci[i]]])){
      
      prot_extractions[[loci[i]]][j]<-strsplit(alignment[[loci[i]]][start[[loci[i]]][j]], " ")
      
      refblock_number[[loci[i]]][j]<-as.numeric(sapply(prot_extractions[[loci[i]]][j], "[", 2))
      
      
      #determines the alignment start by adding -30 to the difference between white spaces found above 
      alignment_start[[loci[i]]]<-refblock_number[[loci[i]]][1]+space_diff[[loci[i]]]
    }
    
    #closes all white space in the alignment file, except for the white space separating the allele and peptide sequence
    alignment[[loci[i]]] <-paste(substr(alignment[[loci[i]]],1,regexpr(" ",text = alignment[[loci[i]]],fixed = TRUE)), gsub(" ","",substr(alignment[[loci[i]]],regexpr(" ",text = alignment[[loci[i]]],fixed = TRUE),nchar(alignment[[loci[i]]]))),sep = "")
    
    #string splits at white spaces to yield allele and peptide sequences
    alignment[[loci[i]]]  <- strsplit(alignment[[loci[i]]]," ", fixed=T)
    
    #binds the previously split strings by row 
    alignment[[loci[i]]] <- do.call(rbind,alignment[[loci[i]]])
    
    #if the pepseq column is equal to the allele column due to premature peptide termination, 
    #insert a blank in place of the allele in the pepseq column 
    alignment[[loci[i]]][which(alignment[[loci[i]]][,1]==alignment[[loci[i]]][,2]),2] <- ""
    
    #renames columns to "alleles" and "pepseq"
    colnames(alignment[[loci[i]]])<-c(paste(loci[[i]], "alleles", sep="_"), "pepseq")
    
    #due to ANHIG formatting, cases where an allele contains newly reference peptide sequences will not 
    #contain the same number of rows as previous reference peptide blocks
    #this for loop is invoked to add "."for all other alleles for each character in the newly reference peptide
    #to preserve structural integrity 
    for(k in 1:length(start[[loci[i]]])){
      if(nrow(alignment[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],])){
        x<-as.data.frame(alignment[[loci[i]]][,1][start[[loci[i]]][1]:end[[loci[i]]][1]][-c(1,2)], stringsAsFactors = F)
        colnames(x)<-paste(loci[[i]], "alleles", sep="_")
        x<-cbind.data.frame(x, pepseq=as.character(paste(rep(".", nchar(tail(alignment[[loci[i]]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
        y<-data.frame(tail(alignment[[loci[i]]],1), stringsAsFactors = F)
        x$pepseq[match(y[,1], x[,1])]<-y$pepseq
        alignment[[loci[i]]]<-as.matrix(rbind(head(alignment[[loci[i]]], -1), x))
        start[[loci[i]]]<-as.numeric(grep("Prot", alignment[[loci[i]]]))
        end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(alignment[[loci[i]]])))}}
    
    #if a locus has extra formatting, resulting in unqeual rows, start and end will be updated to reflect subsetting
    #if a locus has no extra formatting, start and end will remain the same, as procured by earlier code
    for(e in 1:length(start[[loci[i]]])){
      AA_segments[[loci[i]]]<-cbind(AA_segments[[loci[i]]], alignment[[loci[i]]][start[[loci[i]]][e]:end[[loci[i]]][e],])}
    
    #removes first two rows containing AA position and "Prot"
    AA_segments[[loci[i]]] <- AA_segments[[loci[i]]][-c(1,2),]
    
    #designates columns to be combined as every other so allele names are not included
    #in pasting all the amino acid sequences together 
    cols<-seq(0, ncol(AA_segments[[loci[i]]]), by=2)
    AA_segments[[loci[i]]]<-cbind(AA_segments[[loci[i]]][,1], apply(AA_segments[[loci[i]]][,cols], 1 ,paste, collapse = ""))
    
    #creates a new matrix with the number of columns equal to the number of characters in the reference sequence 
    corr_table[[loci[i]]]<-matrix(, nrow = 2, ncol = as.numeric(nchar(AA_segments[[loci[i]]][,2][1])))
    
    #determines alignment length based on the total number of characters plus the alignment start (which is negative ) 
    alignment_length[[loci[i]]]<-as.numeric(nchar(AA_segments[[loci[i]]][,2][1]))+alignment_start[[loci[[i]]]]
    
    #pastes alignment_start to alignment_length together in sequential order, with inDels accounted for 
    #captures output as "w"
    w[[i]] <- capture.output(cat(alignment_start[[loci[i]]]:alignment_length[[loci[i]]]))
    
    #splits string formed by cat for separate character variables
    alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(w[[loci[i]]], " ")))
    
    #eliminates "0", as the alignment sequence from ANHIG does not contain 0
    alignment_positions[[loci[i]]]<-alignment_positions[[loci[i]]][-which(alignment_positions[[loci[i]]] == 0)]
    
    #contains alignment sequence information 
    corr_table[[loci[i]]][2,]<-alignment_positions[[loci[i]]]
    
    #string splits to extract locus in the allele name
    #assigns to new variable "AA_aligned"
    AA_aligned[[loci[i]]]<- as.matrix(do.call(rbind,strsplit(AA_segments[[loci[i]]][,1],"[*]")))
    
    #adds a new column of pasted locus and trimmed two field alleles to AA_aligned
    AA_aligned[[loci[i]]]<- cbind(AA_aligned[[loci[i]]], paste(AA_aligned[[loci[i]]][,1], apply(AA_aligned[[loci[i]]],MARGIN=c(1,2),FUN=GetField,Res=2)[,2], sep="*"))
    
    #binds AA_aligned and AA_segments -- renames columns 
    AA_segments[[loci[i]]] <- cbind(AA_aligned[[loci[i]]], AA_segments[[loci[i]]])
    colnames(AA_segments[[loci[i]]]) <- c("locus", "full_allele", "trimmed_allele", "allele_name", "AAsequence")
    
    #sets refexon to a reference peptide for each HLA locus based on the reference sequences in AA_segments 
    refexon[[loci[i]]] <- rbind(AA_segments[[loci[i]]][1,])[which(rbind(AA_segments[[loci[i]]][1,])[,"locus"]==loci[[i]]),'AAsequence']
    
    #splits AA_sequence column at every amino acid, resulting in a list of split amino acids for each row
    pepsplit[[loci[i]]] <- sapply(AA_segments[[loci[i]]][,"AAsequence"],strsplit,split="*")
    
    #fills in space with NA for alleles with premature termination to make it the same number of characters
    #as the reference sequence 
    pepsplit[[loci[i]]]<- lapply(pepsplit[[loci[i]]],function(x) c(x,rep("NA",nchar(refexon[[loci[i]]])-length(x))))
    
    #binds pep_split together by element in its previous list form by row
    pepsplit[[loci[i]]]<- do.call(rbind,pepsplit[[loci[i]]])
    
    #nullifies row names 
    rownames(pepsplit[[loci[i]]]) <- NULL
    
    #binds all columns together to form desired ouput, as described above
    AA_segments[[loci[i]]] <- cbind.data.frame(AA_segments[[loci[i]]][,1:4],pepsplit[[loci[i]]], stringsAsFactors=FALSE)
    
    #finds positions in AA_segments that have ".", indicating an inDel 
    inDels[[loci[[i]]]]<-colnames(AA_segments[[loci[[i]]]][1, 5:ncol(AA_segments[[loci[[i]]]])][AA_segments[[loci[[i]]]][1, 5:ncol(AA_segments[[loci[[i]]]])] %in% "."])
    
    #inputs AA_segments alignment sequence into the corr_table with "InDel" still present
    corr_table[[loci[[i]]]][1,]<-names(AA_segments[[loci[[i]]]][5:ncol(AA_segments[[loci[[i]]]])])
    
    for(b in 1:length(inDels[[loci[[i]]]])){
      corr_table[[loci[[i]]]][2,][inDels[[loci[[i]]]][[b]]==corr_table[[loci[[i]]]][1,]]<-paste("InDel", b, sep="_")
    }
    
    #fixes enumerations following "InDel"
    corr_table[[loci[[i]]]][2,][!grepl("InDel", corr_table[[loci[[i]]]][2,])]<-(alignment_start[[loci[[i]]]]:((length(corr_table[[loci[[i]]]][2,])-length(corr_table[[loci[[i]]]][2,][grepl("InDel", corr_table[[loci[[i]]]][2,])]))+alignment_start[[loci[[i]]]]))[!(alignment_start[[loci[[i]]]]:((length(corr_table[[loci[[i]]]][2,])-length(corr_table[[loci[[i]]]][2,][grepl("InDel", corr_table[[loci[[i]]]][2,])]))+alignment_start[[loci[[i]]]]))==0]
    
    #renames columns in AA_segments
    colnames(AA_segments[[loci[i]]]) <- c("locus","allele","trimmed_allele","allele_name", 1:ncol(corr_table[[loci[[i]]]]))
    
    #distributes  reference sequence from row 1
    #into all other rows, if they contain a "-"
    #amino acids with changes will not be impacted
    for(k in 5:ncol(AA_segments[[loci[i]]])) {
      AA_segments[[loci[i]]][,k][which(AA_segments[[loci[i]]][,k]=="-")] <- AA_segments[[loci[i]]][,k][1]}  
    
    ##binds a new row to AA_segments at the top of the dataframe -- fills with NAs
    AA_segments[[loci[[i]]]]<-rbind(rep(NA, ncol(AA_segments[[loci[[i]]]])), AA_segments[[loci[[i]]]])
    
    #fills in corresponding alignment sequence in new row
    AA_segments[[loci[[i]]]][1,5:ncol(AA_segments[[loci[[i]]]])]<-corr_table[[loci[[i]]]][2,]
  }
return(AA_segments)
}

##example of AA_segments_maker for HLA
AA_segments<-AA_segments_maker(c("A", "B", "C","DRB1", "DQB1", "DPB1"))


#####motif_finder
#narrows down AA_segments to only alleles that have a given amino acid motif 
motif_finder<-function(motif){

  alignment_corr<-NULL

  #extract loci information 
  loci<-strsplit(motif, "\\*")[[1]][1]
  
  #if AA_segments does not exist (i.e not previously already downloaded and in the local
  #environment) then generate AA_segments)
  if(!exists("AA_segments")){
  #obtains AA_segments df
  AA_segments<-AA_segments_maker(loci)}
  
  #since "DRB" is used as the search criteria for the alignment (IMGTHLA/ANHIG groups all DRB loci
  #into one alignment, AA_segments consists of all DRB loci, not just DRB1)
  #if the loci is DRB1, this conditional statement subsets AA_segments to only DRB1 loci,
  #and if "NA" is present in the locus column for the alignment sequence coordinate row 
  if(loci=="DRB1"){
    AA_segments[[loci]]<-subset(AA_segments[[loci]], (loci==AA_segments[[loci]]$locus) | (is.na(AA_segments[[loci]]$locus)))
  }
  
  #examines motifs to make sure amino acid positions are in the correct order -- sorts numerically
  #if they are not 
  motif<-paste(loci, sep="*",paste(mixedsort(strsplit(strsplit(motif, "*", fixed=T)[[1]][[2]], "~")[[1]]), collapse="~"))
  
  #sets the alignmnet coordinates in the first row of AA_segments to another variable
  alignment_corr[[loci]]<-AA_segments[[loci]][1,]
  
  #for loop for searching amino acid motifs 
  #three total loops are run, where subsequent position*motifs are evaluated against which alleles are present
  #after the previous motif subset 
  #each AA_segments is bound with the alignment coordinates except for the last AA_segments
  for(t in 1:length(strsplit(strsplit(motif, "*", fixed=T)[[1]][[2]], "~")[[1]])){
    AA_segments[[loci]]<-AA_segments[[loci]][which((AA_segments[[loci]][5:ncol(AA_segments[[loci]])][which((str_extract(strsplit(strsplit(motif,"*",fixed=TRUE)[[1]][2],"~",fixed=TRUE)[[1]], "[0-9]+")[[t]]==AA_segments[[loci]][1,5:ncol(AA_segments[[loci]])])==TRUE)]==str_extract(strsplit(strsplit(motif,"*",fixed=TRUE)[[1]][2],"~",fixed=TRUE)[[1]],"[A-Z]")[[t]])==TRUE),]
    
    AA_segments[[loci]]<-rbind(alignment_corr[[loci]], AA_segments[[loci]])}

  
  #if no motifs are found, a warning message is thrown 
  if((nrow(AA_segments[[loci]])==1)){
    warning("Error - zero alleles match this motif. Please try again.")
  }
  
  #if motifs are found, AA_segments[[loci[[i]]]] is returned 
  if((nrow(AA_segments[[loci]])!=1)){
    return(AA_segments[[loci]])}
    
}


#example with actual motif 
X<-motif_finder("DRB1*26F~28E~30Y")

X[c(1,5,9,10), 1:65]

#example with non-existent motif 
motif_finder("DRB1*26F~28E~30Z")


#a function to manipulate the Solberg dataset -- will be adjusted to work with AFND data
dataSubset<-function(dataset, motif){
#reads in Solberg DS
solberg_DS<-as.data.frame(read.delim(dataset), stringsAsFactors=F)

#makes a new column with locus and trimmed allele pasted together named locus_allele
solberg_DS$locus_allele<-paste(solberg_DS$locus, solberg_DS$allele_v3, sep="*")

#orders solberg-DS by population 
solberg_DS<-solberg_DS[order(solberg_DS$popname),]

solberg_DS[,]<-sapply(solberg_DS[, ], as.character)

#subsets the Solberg_DS to only the locus of interest 
solberg_DS<-subset(solberg_DS, solberg_DS$locus==strsplit(motif, "\\*")[[1]][1])

return(solberg_DS)}


#coordinate_converter function
#converts latitude and longtidue to appropriate enumerations based on cardinal directions 
#NORTH and EAST are positive
#N = above equator
#E = east of the prime meridian 
#SOUTH AND WEST ARE NEGATIVE
#S = below equator
#W = west of the prime meridian 
coordinate_converter<-function(heatmapdata){
  #latitude conversions
  if(any((grepl("S", heatmapdata$latit))==TRUE)){
    heatmapdata$latit[which((grepl("S", heatmapdata$latit))==TRUE)]<-as.numeric(paste("-", str_extract(heatmapdata$latit[which((grepl("S", heatmapdata$latit))==TRUE)],"\\-*\\d+\\.*\\d*"), sep=""))}
  
  if(any((grepl("N", heatmapdata$latit))==TRUE)){
    heatmapdata$latit[which((grepl("N", heatmapdata$latit))==TRUE)]<-as.numeric(str_extract(heatmapdata$latit[which((grepl("N", heatmapdata$latit))==TRUE)],"\\-*\\d+\\.*\\d*"))
  }
  
  #longitude conversions
  if(any((grepl("W", heatmapdata$longit))==TRUE)){
    heatmapdata$longit[which((grepl("W", heatmapdata$longit))==TRUE)]<-as.numeric(paste("-", str_extract(heatmapdata$longit[which((grepl("W", heatmapdata$longit))==TRUE)],"\\-*\\d+\\.*\\d*"), sep=""))}
  
  if(any((grepl("E", heatmapdata$longit))==TRUE)){
    heatmapdata$longit[which((grepl("E", heatmapdata$longit))==TRUE)]<-as.numeric(str_extract(heatmapdata$longit[which((grepl("E", heatmapdata$longit))==TRUE)],"\\-*\\d+\\.*\\d*"))
  }
  return(heatmapdata)}


#PALM () function  -- Population Allele Locating Mapmaker 
##default for color is set to TRUE,
#where color defines whether the user wants grayscale or colorscale heatmap 
##default for filter_migrant is set to TRUE, where migrant admixed populations and OTH populations
#are filtered out
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

#example of PALM() where output has color
PALM("1-locus-alleles.dat", "DRB1*26F~28E~30Y")

#example of PALM() where output has no color
PALM("1-locus-alleles.dat", "DRB1*26F~28E~30Y", color=FALSE)


