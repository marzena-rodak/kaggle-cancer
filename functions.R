
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

catch_mutation_grams <- function(sentence){
  x <- regmatches(sentence,gregexpr('\\w+\\s(mutation)',sentence))
  return(x)
}

sentence_with_variation <- function(sent_list,var){
  
  sent_var <- list()
  
  for(k in 1:length(sent_list)){
    if(grepl(var,sent_list[[k]])){
      sent_var <- c(sent_var,sent_list[[k]])
    }
  }
  
  return(sent_var)
}

sentence_with_mutat <- function(sent_list){
  
  sent_mut <- list()
  
  for(k in 1:length(sent_list)){
    if(grepl('\\s(mutat)',sent_list[[k]])){
      sent_mut <- c(sent_mut,sent_list[[k]])
    }
  }
  
  return(sent_mut)
}