# Load CSV files
cat("Read data")
train_text <- do.call(rbind,strsplit(readLines('Data/training_text'),'||',fixed=T))
train_text <- as.data.table(train_text)
train_text <- train_text[-1,]
colnames(train_text) <- c("ID", "Text")
train_text$ID <- as.numeric(train_text$ID)

test_text <- do.call(rbind,strsplit(readLines('Data/test_text'),'||',fixed=T))
test_text <- as.data.table(test_text)
test_text <- test_text[-1,]
colnames(test_text) <- c("ID", "Text")
test_text$ID <- as.numeric(test_text$ID)

train <- fread("../input/training_variants", sep=",", stringsAsFactors = T)
test <- fread("../input/test_variants", sep=",", stringsAsFactors = T)
train <- merge(train,train_text,by="ID")
test <- merge(test,test_text,by="ID")
rm(test_text,train_text);gc()

test$Class <- -1
data <- rbind(train,test)
rm(train,test);gc()
