#'Converts longitude and latitude with cardinal directions to proper coordinates
#'
#'The Solberg dataset contains some coordinates associated with a cardinal direction, such as 50S. This function converts those coordinates to their enumerations -- in this case, 50S is -50. North (above equator) and East (east of prime meridian) positions are positive, and South (below equator) and West (west of prime meridian) are negative.
#'
#'@importFrom stringr str_extract
#'
#'@param heatmapdata Manipulated Solberg dataset produced from dataSubset.
#'
#'@return Solberg dataset with appropriate latitude and longitude enumerations
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
