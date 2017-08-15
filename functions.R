
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
  # x <- gsub('[A-z0-9]+[)]','',x)
  # x <- gsub('[A-z0-9]+[]]','',x)
  # x <- gsub('[A-z0-9]+[}]','',x)
  return(x)
}

preprocess_Corpus <-function(txt){
  txt <- Corpus(VectorSource(txt))
  txt <- tm_map(txt, stripWhitespace)
  txt <- tm_map(txt, content_transformer(tolower))
  txt <- tm_map(txt, removeReferences)
  txt <- tm_map(txt, removeReferences)
  txt <- tm_map(txt, removeReferences)
  txt <- tm_map(txt, removeNumbers2)
  txt <- tm_map(txt, removePunctuation, preserve_intra_word_dashes = TRUE)
  txt <- tm_map(txt, removeWords, c(stopwords("english"),'figure','fig','table','plot','chart','graph'))
  txt <- tm_map(txt, stripWhitespace)
  txt <- tm_map(txt, stemDocument, language="english")
  return(txt)
}

catch_mutation_grams <- function(sentence){
  x <- regmatches(sentence,gregexpr('\\w+\\s(mutation)',sentence))
  return(x)
}

sentence_with_keyword <- function(sent_list,keyword){
  
  sent_keyword <- list()
  
  if( !is.na(keyword) ){
    for(k in 1:length(sent_list)){
      if(grepl(paste0("(",paste(keyword,collapse="|"),")"), sent_list[[k]], ignore.case = TRUE)){
        sent_keyword <- c(sent_keyword,sent_list[[k]])
      }
    }
  }
  
  return(sent_keyword)
}


convert_text_to_sentences <- function(text, lang = "en") {
  # Ensure there are spaces after each dot.
  text <- stripWhitespace(gsub('\\.', '. ', text))
  
  # Function to compute sentence annotations using the Apache OpenNLP Maxent sentence detector employing the default model for language 'en'. 
  sentence_token_annotator <- Maxent_Sent_Token_Annotator(language = lang)
  
  # Convert text to class String from package NLP
  text <- as.String(text)
  
  # Sentence boundaries in text
  sentence.boundaries <- annotate(text, sentence_token_annotator)
  
  # Extract sentences
  sentences <- text[sentence.boundaries]
  
  # return sentences
  return(sentences)
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

LogLossSummary <- function (data, lev = NULL, model = NULL) {
  LogLos <- function(actual, pred, eps = 1e-15) {
    stopifnot(all(dim(actual) == dim(pred)))
    pred[pred < eps] <- eps
    pred[pred > 1 - eps] <- 1 - eps
    -sum(actual * log(pred)) / nrow(pred) 
  }
  if (is.character(data$obs)) data$obs <- factor(data$obs, levels = lev)
  pred <- data[, "pred"]
  obs <- data[, "obs"]
  isNA <- is.na(pred)
  pred <- pred[!isNA]
  obs <- obs[!isNA]
  data <- data[!isNA, ]
  cls <- levels(obs)
  if (length(obs) + length(pred) == 0) {
    out <- rep(NA, 2)
  } else {
    pred <- factor(pred, levels = levels(obs))
    require("e1071")
    out <- unlist(e1071::classAgreement(table(obs, pred)))[c("diag", "kappa")]
    probs <- data[, cls]
    actual <- model.matrix(~ obs - 1)
    out2 <- LogLos(actual = actual, pred = probs)
  }
  out <- c(out, out2)
  names(out) <- c("Accuracy", "Kappa", "LogLoss")
  if (any(is.nan(out))) out[is.nan(out)] <- NA
  out
}


processCorpus <- function (corpus) {
  corpus <- tm_map(corpus, stripWhitespace)
  corpus <- tm_map(corpus, content_transformer(tolower))
  corpus <- tm_map(corpus, stemDocument, language="english")
  corpus <- tm_map(corpus, removePunctuation, preserve_intra_word_dashes = TRUE)
  corpus <- tm_map(corpus, removeWords, stopwords("english")) 
  corpus <- tm_map(corpus, function (x) {
    gsub("\\s*(?<!\\B|-)\\d+(?!\\B|-)\\s*", "", x, perl = TRUE)})
  corpus <- tm_map(corpus, PlainTextDocument)
  return (corpus)
}

getDictionary <- function (filePath) {
  if (is.null(filePath)) return(NULL)
  ret <- do.call(rbind, strsplit(readLines(filePath), "\t", fixed = TRUE))
  ret <- setNames(data.table(ret, stringsAsFactors = FALSE), 
                  c("code", "concept_name", "parents", 
                    "synonyms", "definition", "display_name",
                    "concept_status", "semantic_type"))
  
  corpus <- Corpus(VectorSource(gsub("\\|", " ", ret$synonyms)))
  corpus <- processCorpus (corpus)
  
  terms <- DocumentTermMatrix(
    corpus, 
    control = list(
      minWordLength = 3,
      weighting = function(x) weightTfIdf(x, normalize = FALSE)))
  
  return (as.vector(terms$dimnames$Terms))
}