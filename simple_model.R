train_cont <- create_container(dtm, labels = train$Class, trainSize = 1:300, virgin = FALSE)
train_mat <- as.matrix(train_cont@training_matrix)
colnames(train_mat) <- dtm$dimnames[[2]]
train_lab <- as.factor(paste0("c", train$Class[1:300]))

trControl <- trainControl(method = "cv", number = 3, classProbs = TRUE, summaryFunction = LogLossSummary)
xgb.grid <- expand.grid(nrounds = c(50, 100), max_depth = c(5, 7), eta = 0.1, gamma = 1,
                        colsample_bytree = 0.1, min_child_weight = 2, subsample = 1)
system.time(
  caret.fit <- train(train_mat, train_lab, method = 'xgbTree', metric = 'LogLoss', trControl = trControl, tuneGrid = xgb.grid, maximize = FALSE)
)
caret.fit
