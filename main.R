# Install and load libraries ----------------------------------------------
# install.packages("openNLPmodels.en", repos = "http://datacube.wu.ac.at/", type = "source")
library(caret) ### IMPORTANT! Load caret before NLP; `annotate` from NLP is masked.
library(NLP)
library(openNLP)
library(tm)
library(RTextTools)
library(RWeka)
library(data.table)
library(SnowballC)
library(xgboost)
library(stringdist)
library(dplyr)
library(slam)
source('functions.R')
library(foreach)
library(doParallel)

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

txt <- preprocess_Corpus(data$Text)
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
sums <- colapply_simple_triplet_matrix(dtm,FUN=sum)
sums <- sort(sums, decreasing=T)
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



# Variation distance analysis ---------------------------------------------



# sum_var <- train %>% group_by(Variation) %>% summarize(count = n()) %>% filter(count > 2) %>% arrange(-count)
sum_var <- train %>% group_by(Variation) %>% summarize(count = n()) %>% arrange(-nchar(as.character(Variation)))
x <- train[train[,'Variation'] %in% sum_var$Variation, c('Variation','Class')] %>% group_by(Variation, Class) %>% summarize(count = n())
ggplot(x) + geom_point( aes(x=Variation, y=Class, size = count))

# Get word and n-gram frequency for Variations
var_corp <- VCorpus(VectorSource(as.character(train$Variation)))
var_corp <- tm_map(var_corp, content_transformer(tolower))
var_corp <- tm_map(var_corp, stemDocument, language="english")

# Singles
var_corp_1g <- DocumentTermMatrix(var_corp)

# Bigrams
options(mc.cores=1)
var_corp_2g <- DocumentTermMatrix(var_corp, control=list(tokenize=BigramTokenizer, wordLengths = c(1, Inf)))

# To get the bigram dist, we use the slam package for ops with simple triplet mat
sums.1g <- colapply_simple_triplet_matrix(var_corp_1g,FUN=sum)
sums.1g <- sort(sums.1g, decreasing=T)

# To get the bigram dist, we use the slam package for ops with simple triplet mat
sums.2g <- colapply_simple_triplet_matrix(var_corp_2g,FUN=sum)
sums.2g <- sort(sums.2g, decreasing=T)

df_plot.1g <- data.frame(n_gram = factor(names(sums.1g), ordered = TRUE, levels = names(sums.1g)), freq = sums.1g)
df_plot.2g <- data.frame(n_gram = factor(names(sums.2g), ordered = TRUE, levels = names(sums.2g)), freq = sums.2g)
gg.1g <- ggplot(df_plot.1g[1:20, ], aes(x = n_gram, y = freq)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg.2g <- ggplot(df_plot.2g[1:20, ], aes(x = n_gram, y = freq)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# gg.1g
# gg.2g

## Grep different mutations and transform them to one-hot encoding

# Identify and encode deletions
train$is_del <- ifelse(grepl("del", train$Variation, ignore.case = T), 1, 0)

# Identify and encode insertions
train$is_ins <- ifelse(grepl("ins", train$Variation, ignore.case = T), 1, 0)

# Identify and encode fusions
train$is_fus <- ifelse(grepl("fus", train$Variation, ignore.case = T), 1, 0)

# Identify and encode truncation
train$is_trunc <- ifelse(grepl("trunc", train$Variation, ignore.case = T), 1, 0)

# Identify and encode methylations
train$is_methyl <- ifelse(grepl("methyl", train$Variation, ignore.case = T), 1, 0)

# Identify and encode amplifications
train$is_amp <- ifelse(grepl("amp", train$Variation, ignore.case = T), 1, 0)

# Identify and encode silencing
train$is_sil <- ifelse(grepl("sil", train$Variation, ignore.case = T), 1, 0)

# Identify and encode overexpression
train$is_expr <- ifelse(grepl("expr", train$Variation, ignore.case = T), 1, 0)

# Identify and encode splicing
train$is_splice <- ifelse(grepl("splice", train$Variation, ignore.case = T), 1, 0)

# Identify and encode exon variations
train$is_exon <- ifelse(grepl("exon", train$Variation, ignore.case = T), 1, 0)





######Building dtm separetly from smaller number of sentences

#sub_train_full <- train[!is.na(train$First_Var),]

index <- createDataPartition(train$Class, times = 1, p=0.01, list = FALSE)
sub_train <- train[index,]


# Gene variation extraction -----------------------------------------------
subtrain_Var<-list()

mut_vec<-c("delet","insert",'fusion',"truncat","methyl","amplif","silenc","express","splic","exon")
pom_Var<- apply(sub_train[,c("is_del","is_ins","is_fus","is_trunc","is_methyl","is_amp","is_sil","is_expr","is_splice","is_exon")],
                1, function(x) mut_vec[x==1])

r <-gregexpr('[A-Z0-9]{2,}', sub_train$Variation)

for(j in 1:nrow(sub_train)){
  # if(r[[j]][1]!=-1){
    subtrain_Var[[j]] <- c(regmatches(sub_train$Variation[j],r[j])[[1]], pom_Var[[j]])
  # }
}

cl <- makeCluster(4)
registerDoParallel(cl)

gene_doc_list = foreach(p = 1:nrow(sub_train),.packages = c("NLP","openNLP","tm")) %dopar% {
  print(p)
  text <- sub_train$Text[p]
  sent <- convert_text_to_sentences(text)
  sent_gene <- sentence_with_keyword(sent,sub_train$Gene[p])
  paste(sent_gene, collapse=" ")
}

var_doc_list = foreach(p = 1:nrow(sub_train),.packages = c("NLP","openNLP","tm")) %dopar% {
  print(p)
  text <- sub_train$Text[p]
  sent <- convert_text_to_sentences(text)
  sent_var <- sentence_with_keyword(sent, subtrain_Var[[p]])
  paste(sent_var, collapse=" ")
}

stopCluster(cl)

gene_txt <- preprocess_Corpus(gene_doc_list)
gene_txt <- tm_map(gene_txt, PlainTextDocument)
var_txt <- preprocess_Corpus(var_doc_list)
var_txt <- tm_map(var_txt, PlainTextDocument)

#One-gram term matrix
gene_dtm <- DocumentTermMatrix(gene_txt, control = list(weighting = weightTfIdf))
gene_dtm <- removeSparseTerms(gene_dtm, 0.95)

var_dtm <- DocumentTermMatrix(var_txt, control = list(weighting = weightTfIdf))
var_dtm <- removeSparseTerms(var_dtm, 0.95)

#Bigram term matrix
gene_dtm_2 <- DocumentTermMatrix(gene_txt, control=list(tokenize=BigramTokenizer, wordLengths = c(1, Inf)))
gene_dtm_2 <- removeSparseTerms(gene_dtm_2, 0.95)

var_dtm_2 <- DocumentTermMatrix(var_txt, control=list(tokenize=BigramTokenizer, wordLengths = c(1, Inf)))
var_dtm_2 <- removeSparseTerms(var_dtm_2, 0.95)

gene_dtm_comb <- c(gene_dtm,gene_dtm_2)
var_dtm_comb <- c(var_dtm,var_dtm_2)

#######Build even smaller dtm using NCLt dictionary

#One-gram term matrix
nci_gene_dtm <- DocumentTermMatrix(gene_txt, control = list(dictionary = medDico, weighting = weightTfIdf))
nci_gene_dtm <- removeSparseTerms(nci_gene_dtm, 0.95)

nci_var_dtm <- DocumentTermMatrix(var_txt, control = list(dictionary = medDico, weighting = weightTfIdf))
nci_var_dtm <- removeSparseTerms(nci_var_dtm, 0.95)



# NCI Thesaurus -----------------------------------------------------------

medDico <- getDictionary("Data/Thesaurus.txt")

text <- sub_train$Text[empty_gene[1]]
key <- sub_train$Gene[empty_gene[1]]
sub_train$ID[empty_gene[1]]
sub_train$ID[empty_var]
