##PAL (Population Allele Locater)
#v 0.3
#By: Liv Tran
#5/21/19


###old script to get AA_segments
##REQUIRED PACKAGES 
require(data.table)
require(stringr)
require(gtools)

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

loci="DRB1"


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


##example of AA_segments_maker for HL
AA_segments<-AA_segments_maker(c("A", "B", "C","DRB1"))

#### start script for motif_finder
motif_finder<-function(motif){

  #extract loci information 
  loci<-strsplit(motif, "\\*")[[1]][1]
  
  AA_segments<-AA_segments_maker(loci)
  
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
    
    if(t!=3){
    AA_segments[[loci]]<-rbind(alignment_corr[[loci]], AA_segments[[loci]])}
  }
  
  
  
  #if no motifs are found, a warning message is thrown 
  if((nrow(AA_segments[[loci]])==0)){
    warning("Error - zero alleles match this motif. Please try again.")
  }
  
  #if motifs are found, AA_segments[[loci[[i]]]] is returned 
  if((nrow(AA_segments[[loci]])!=0)){
    return(AA_segments[[loci]])}
    
}

####EXAMPLES


#example with actual motif 
motif_finder("DRB1*26F~28E~30Y")

#example with non-exisent motif 
motif_finder("DRB1*26F~28E~30Z")



#Population Allele Locator (PAL) function
#An independent function for getting allele frequencies from a population for a given motif 
#for the Solberg dataset 
PAL<-function(dataset, motif){
  
  #saved to Alleles With Motifs (AWMs)
  AWMs<-motif_finder(motif)
  
  #reads in Solberg DS
  solberg_DS<-as.data.frame(read.delim(dataset), stringsAsFactors=F)
  
  #makes a new column with locus and trimmed allele pasted together named locus_allele
  solberg_DS$locus_allele<-paste(solberg_DS$locus, solberg_DS$allele_v3, sep="*")
  
  #orders solberg-DS by population 
  solberg_DS<-solberg_DS[order(solberg_DS$popname),]
  
  solberg_DS[,]<-sapply(solberg_DS[, ], as.character)
  
  #makes an empty list named unique_AWMs, where the name of each element is after a unique AWM
  unique_AWMs<-sapply(unique(AWMs$trimmed_allele), function(x) NULL)
  
  #finds unique_AWMs in Solberg dataset and extracts the allele frequency and locus_allele column  
  for(y in 1:length(unique_AWMs)){
    unique_AWMs[[y]]<-solberg_DS[,c(2,3,5,6,11,14)][solberg_DS$locus_allele %in% names(unique_AWMs[y]),]}
  
  #subsets out locus_allele pairs with the motif from the alignment but aren't present in the Solberg ds
  unique_AWMs<-unique_AWMs[sapply(unique_AWMs, nrow)>0]
  
  #sets population names to a variable named popnames to NULL 
  popnames<-NULL
  
  #sets "to be mapped dataset" to tbm_ds
  tbm_ds<-NULL
  
  #inserts all popnames column from unique_AWMs to popnames
  #inserts popnames, continent, latitude, and longitude information into tbm_ds
  for(i in 1:length(unique_AWMs)){
    popnames[[i]]<-unique_AWMs[[i]][1]
    tbm_ds[[i]]<-unique_AWMs[[i]][c(1,2,3,4)]}
  
  #unlists popnames -- finds unique popnames 
  popnames<-unique(unlist(popnames))
  
  #melts tbm_ds data into one data frame
  tbm_ds<-melt(tbm_ds, id.vars=c("popname", "contin", "latit", "longit"))
  
  #removes duplicates
  tbm_ds<-tbm_ds[,c(1:ncol(tbm_ds))][!duplicated(tbm_ds$popname),]
  
  #removes L1 column based on melt 
  tbm_ds$L1<-NULL
  
  #creates a variable named Population Allele Frequencies (PAF), where each element is named after a 
  #unique popname
  #makes each element a list to take in allele frequencies in the next for loop 
  PAF<-sapply(popnames, function(x) list())
  
  #for loop for finding allele frequencies from unique_AWMs if a PAF name matches a unique_AWM element
  for(i in 1:length(PAF)){
    for(j in 1:length(unique_AWMs)){
      
      #matches each PAF name to each unique_AWMs element, extracts allele frequencies from each element
      #if that name is not found for a unique_AWMs element, "NA" is present in place of a value
      PAF[[i]][j]<-unique_AWMs[[j]][,5][match(names(PAF[i]), unique_AWMs[[j]][,1])]
    }
    
    #unlists each PAF element, and subsets out any NAs
    #if a PAF name has duplicates, the allele frequencies are summed
    #PAF names with only one allele value are used and unaffected by the sum() function 
    PAF[[i]]<-sum(as.numeric(unlist(PAF[[i]])[!is.na(unlist(PAF[[i]]))]))
  }
  
  #melts PAF into a two columned df
  PAF<-melt(PAF)
  
  #renames column names 
  colnames(PAF)<-c("allele_freq", "popname")
  
  #merges tbm_ds with PAF information
  tbm_ds<-merge(tbm_ds, PAF, by="popname")
  
  return(tbm_ds)
}
  

#example of PAL()
#saved to Heat Map Data (HMD)
HMD<-PAL("1-locus-alleles.dat", "DRB1*26F~28E~30Y") 

