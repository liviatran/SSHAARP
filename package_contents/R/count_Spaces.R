#'Counts whitespaces in a string
#'
#'Provided string is split, and white spaces are counted in between split elements. Written by Josh Bredeweg (reddit user jbraids1421).
#'@param x A string.
#'
#'@return Number of white spaces between elements of a split string.
#'
#'@export
#'
#'@examples countSpaces("The quick  brown   fox jumped    over  the lazy   dog.")
#
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
}
