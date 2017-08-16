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
library(plyr)
library(dplyr)
library(slam)
library(foreach)
library(doParallel)
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

# Grep different mutations and transform them to one-hot encoding ---------

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

index <- createDataPartition(train$Class, times = 1, p=0.4, list = FALSE)
sub_train <- train[index,]
# NCI Thesaurus
medDico <- getDictionary("Data/Thesaurus.txt")

# Gene variation extraction -----------------------------------------------
subtrain_Var<-list()

mut_vec<-c("delet","insert",'fusion',"truncat","methyl","amplif","silenc","express","splic","exon")
mut_col <- c("is_del","is_ins","is_fus","is_trunc","is_methyl","is_amp","is_sil","is_expr","is_splice","is_exon")
pom_Var<- apply(sub_train[,mut_col],
                1, function(x) mut_vec[x==1])

r <-gregexpr('[A-Z0-9]{2,}', sub_train$Variation)

for(j in 1:nrow(sub_train)){
  # if(r[[j]][1]!=-1){
    subtrain_Var[[j]] <- c(regmatches(sub_train$Variation[j],r[j])[[1]], pom_Var[[j]])
  # }
}

zero_ind <- sapply(subtrain_Var, length) == 0

# Remove observations with no Variation
sub_train <- sub_train[!zero_ind, ]
subtrain_Var <- subtrain_Var[!zero_ind]

cl <- makeCluster(4)
registerDoParallel(cl)

gene_doc_list = foreach(p = 1:nrow(sub_train),.packages = c("NLP","openNLP","tm")) %dopar% {
  text <- sub_train$Text[p]
  sent <- convert_text_to_sentences(text)
  sent_gene <- sentence_with_keyword(sent,sub_train$Gene[p])
  paste(sent_gene, collapse=" ")
}

var_doc_list = foreach(p = 1:nrow(sub_train),.packages = c("NLP","openNLP","tm")) %dopar% {
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

gene_dtm_comb <- cbind(gene_dtm, gene_dtm_2, as.matrix(sub_train[, mut_col]))
gene_dtm_comb$dimnames[[2]] <- c(gene_dtm$dimnames[[2]], gene_dtm_2$dimnames[[2]], mut_col)
gene_dtm_ind <- apply(gene_dtm_comb , 1, sum) > 0
gene_dtm_comb <- gene_dtm_comb[gene_dtm_ind, ]
gene_dtm_comb_lab <- sub_train$Class[gene_dtm_ind]

var_dtm_comb <- cbind(var_dtm, var_dtm_2, as.matrix(sub_train[, mut_col]))
var_dtm_comb$dimnames[[2]] <- c(var_dtm$dimnames[[2]], var_dtm_2$dimnames[[2]], mut_col)
var_dtm_ind <- apply(var_dtm_comb , 1, sum) > 0
var_dtm_comb <- var_dtm_comb[var_dtm_ind, ]
var_dtm_comb_lab <- sub_train$Class[var_dtm_ind]

#######Build even smaller dtm using NCLt dictionary

#One-gram term matrix
nci_gene_dtm <- DocumentTermMatrix(gene_txt, control = list(dictionary = medDico, weighting = weightTfIdf))
nci_gene_dtm <- removeSparseTerms(nci_gene_dtm, 0.95)
nci_gene_dtm <- cbind(nci_gene_dtm, as.matrix(sub_train[, mut_col]))
nci_gene_ind <- apply(nci_gene_dtm , 1, sum) > 0
nci_gene_dtm <- nci_gene_dtm[nci_gene_ind, ]
nci_gene_dtm_lab <- sub_train$Class[nci_gene_ind]

nci_var_dtm <- DocumentTermMatrix(var_txt, control = list(dictionary = medDico, weighting = weightTfIdf))
nci_var_dtm <- removeSparseTerms(nci_var_dtm, 0.95)
nci_var_dtm <- cbind(nci_var_dtm, as.matrix(sub_train[, mut_col]))
nci_var_ind <- apply(nci_var_dtm , 1, sum) > 0
nci_var_dtm <- nci_var_dtm[nci_var_ind, ]
nci_var_dtm_lab <- sub_train$Class[nci_var_ind]

## We have following DTMs:
# gene_dtm, car_dtm, gene_dtm_2, var_dtm_2, gene_dtm_comb, var_dtm_comb, nci_gene_dtm, nci_var_dtm
## We will use four of them
# gene_dtm_comb, var_dtm_comb, nci_gene_dtm, nci_var_dtm

# text <- sub_train$Text[empty_gene[1]]
# key <- sub_train$Gene[empty_gene[1]]
# sub_train$ID[empty_gene[1]]
# sub_train$ID[empty_var]
