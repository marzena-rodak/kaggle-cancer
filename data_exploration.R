# summary(train)
# levels(train$Variation)

# Basic TF-IDF for Text ---------------------------------------------------

index <- createDataPartition(train$Class, times = 1, p=0.1, list = FALSE)
data <- train[index,]

txt <- preprocess_Corpus(data$Text)
dtm <- DocumentTermMatrix(txt, control = list(weighting = weightTfIdf))
dtm <- removeSparseTerms(dtm, 0.95)
#data_dtm <- cbind(data, as.matrix(dtm))

#Add n-grams to dtm

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