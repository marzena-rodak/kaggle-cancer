
# User-defined functions --------------------------------------------------

read <- function(fileName) {
  data <- readLines(con <- file(fileName))
  close(con)
  dataframe <- as.data.frame(data,stringsAsFactors = FALSE)
  dataframe <- dataframe[2:nrow(dataframe),,drop=FALSE]
  dataframe$ID <- ""
  
  for(k in 1:nrow(dataframe)){
    pom <- attr(regexpr('\\d+',dataframe[k,1]),'match.length')
    dataframe[k,2] <-substring(dataframe[k,1],1,pom)
    dataframe[k,1] <- substring(dataframe[k,1],pom+3)
  }
  
  colnames(dataframe) <- c('Text','ID')
  dataframe$ID <- as.numeric(dataframe$ID)
  
  return(dataframe[,c(2,1)])
}

removeNumbers2 <-function(x){
  x <- gsub('(^|\\s)\\d+(\\s|\\.|$)',' ',x)
  return(x)
}

removeReferences <-function(x){
  x <- gsub('[(][^()]+[)]','',x) #remove everything in parenthesis ()
  x <- gsub('[[][^()]+[]]','',x) #remove everything in parenthesis []
  x <- gsub('[{][^()]+[}]','',x) #remove everything in parenthesis {}
  x <- gsub('[A-z0-9]+[)]','',x)
  x <- gsub('[A-z0-9]+[]]','',x)
  x <- gsub('[A-z0-9]+[}]','',x)
  return(x)
}