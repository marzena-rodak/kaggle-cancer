# Install and load libraries ----------------------------------------------
library(data.table)
library(tm)
library(SnowballC)
library(caret)
library(xgboost)
library(RTextTools)
library(RWeka)
library(stringdist)
library(dplyr)
source('functions.R')

# Load data ---------------------------------------------------------------

train_var <- read.csv('Data/training_variants')
system.time({
train_text <- read('Data/training_text')
})
test_var <- read.csv('Data/test_variants')
test_text <- read('Data/test_text')
train <- merge(train_var,train_text,by="ID")
test <- merge(test_var,test_text,by="ID")
rm(test_text,train_text,test_var,train_var); gc()

train$Class <- as.factor(train$Class)
summary(train)
levels(train$Variation)

# Basic TF-IDF for Text ---------------------------------------------------

index <- createDataPartition(train$Class, times = 1, p=0.1, list = FALSE)
data <- train[index,]

txt <- Corpus(VectorSource(data$Text))
txt <- tm_map(txt, stripWhitespace)
txt <- tm_map(txt, content_transformer(tolower))
txt <- tm_map(txt, removeReferences)
txt <- tm_map(txt, removeNumbers2)
txt <- tm_map(txt, removePunctuation, preserve_intra_word_dashes = TRUE)
txt <- tm_map(txt, removeWords, c(stopwords("english"),'figure','fig','table','plot','chart','graph'))
txt <- tm_map(txt, stripWhitespace)
txt <- tm_map(txt, stemDocument, language="english")
dtm <- DocumentTermMatrix(txt, control = list(weighting = weightTfIdf))
dtm <- removeSparseTerms(dtm, 0.95)
#data_dtm <- cbind(data, as.matrix(dtm))

#Add n-grams to dtm
BigramTokenizer <- function(x) {RWeka::NGramTokenizer(x, RWeka::Weka_control(min=2, max=2))}
ThreegramTokenizer <- function(x) {RWeka::NGramTokenizer(x, RWeka::Weka_control(min=3, max=3))}
FourgramTokenizer <- function(x) {RWeka::NGramTokenizer(x, RWeka::Weka_control(min=4, max=4))}

# Bigrams
options(mc.cores=1)
dtm.2g <- DocumentTermMatrix(txt, control=list(tokenize=BigramTokenizer,weighting = weightTfIdf))

#Threegrams
options(mc.cores=1)
dtm.3g <- DocumentTermMatrix(txt, control=list(tokenize=ThreegramTokenizer,weighting = weightTfIdf))

#Fourgrams
options(mc.cores=1)
dtm.4g <- DocumentTermMatrix(txt, control=list(tokenize=FourgramTokenizer,weighting = weightTfIdf))




# # To get the word dist, we use the slam package for ops with simple triplet mat
# library(slam)
# sums <- colapply_simple_triplet_matrix(dtm,FUN=sum)
# sums <- sort(sums, decreasing=T)
# 
# # To get the bigram dist, we use the slam package for ops with simple triplet mat
# sums.2g <- colapply_simple_triplet_matrix(dtm.2g,FUN=sum)
# sums.2g <- sort(sums.2g, decreasing=T)
# 
# sums.3g <- colapply_simple_triplet_matrix(dtm.3g,FUN=sum)
# sums.3g <- sort(sums.3g, decreasing=T)
# 
# sums.4g <- colapply_simple_triplet_matrix(dtm.4g,FUN=sum)
# sums.4g <- sort(sums.4g, decreasing=T)



# Variation distance analysis ------------------------------------------------------



sum_var <- train %>% group_by(Variation) %>% summarize(count = n()) %>% filter(count > 2) %>% arrange(-count)
x <- train[train[,'Variation'] %in% sum_var$Variation, c('Variation','Class')] %>% group_by(Variation, Class) %>% summarize(count = n())
ggplot(x) + geom_point( aes(x=Variation, y=Class, size = count))

