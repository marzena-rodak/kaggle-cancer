train_cont <- prep_container(gene_dtm_comb, gene_dtm_comb_lab)

trControl <- trainControl(method = "cv", number = 3, classProbs = TRUE, summaryFunction = LogLossSummary)
xgb.grid <- expand.grid(nrounds = 75,
                        max_depth = c(9, 13, 17),
                        eta = c(0.07, 0.15),
                        gamma = c(0.5, 1),
                        colsample_bytree = 0.1,
                        min_child_weight = c(2),
                        subsample = 0.5)

system.time(
  caret.fit <- train(train_cont$X, train_cont$y, method = 'xgbTree', metric = 'LogLoss',
                     trControl = trControl, tuneGrid = xgb.grid, maximize = FALSE)
)

caret.fit



# var_dtm_comb ------------------------------------------------------------

# eXtreme Gradient Boosting 
# 
# 1083 samples
#  608 predictor
#    9 classes: 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9' 
# 
# No pre-processing
# Resampling: Cross-Validated (3 fold) 
# Summary of sample sizes: 722, 721, 723 
# Resampling results across tuning parameters:
# 
#   max_depth  gamma  min_child_weight  Accuracy   Kappa      LogLoss 
#   5           1     2                 0.5253925  0.3744265  1.410027
#   5           1     4                 0.5180031  0.3640088  1.437010
#   5           1     8                 0.4672206  0.2919919  1.499601
#   5           4     2                 0.4921642  0.3224267  1.494055
#   5           4     4                 0.4940058  0.3229246  1.508756
#   5           4     8                 0.4607621  0.2763053  1.558619
#   5          12     2                 0.4469141  0.2367933  1.624388
#   5          12     4                 0.4090587  0.1780333  1.649953
#   5          12     8                 0.4155222  0.1933581  1.663730
#   7           1     2                 0.5290784  0.3777042  1.389831 (!)
#   7           1     4                 0.5059993  0.3486841  1.436978
#   7           1     8                 0.4773622  0.3051638  1.495402
#   7           4     2                 0.5124706  0.3482056  1.488375
#   7           4     4                 0.4894120  0.3142959  1.505096
#   7           4     8                 0.4496945  0.2607541  1.557500
#   7          12     2                 0.4294009  0.2143306  1.631641
#   7          12     4                 0.4238607  0.2058376  1.650100
#   7          12     8                 0.4053832  0.1755692  1.670897
#   9           1     2                 0.5355521  0.3877620  1.396976
#   9           1     4                 0.5069228  0.3480886  1.426698
#   9           1     8                 0.4783112  0.3064810  1.499835
#   9           4     2                 0.4902996  0.3218624  1.493323
#   9           4     4                 0.4958602  0.3245885  1.497688
#   9           4     8                 0.4570686  0.2706093  1.567204
#   9          12     2                 0.4266077  0.2107666  1.631526
#   9          12     4                 0.4210752  0.2013279  1.643607
#   9          12     8                 0.4136730  0.1911971  1.669206
# 
# Tuning parameter 'nrounds' was held constant at a value of 75
# Tuning parameter 'eta' was held
#  constant at a value of 0.07
# Tuning parameter 'colsample_bytree' was held constant at a value of
#  0.1
# Tuning parameter 'subsample' was held constant at a value of 0.5
# LogLoss was used to select the optimal model using  the smallest value.
# The final values used for the model were nrounds = 75, max_depth = 7, eta = 0.07, gamma =
#  1, colsample_bytree = 0.1, min_child_weight = 2 and subsample = 0.5.

# gene_dtm_comb -----------------------------------------------------------

# eXtreme Gradient Boosting 
# 
# 1242 samples
# 2721 predictors
#    9 classes: 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9' 
# 
# No pre-processing
# Resampling: Cross-Validated (3 fold) 
# Summary of sample sizes: 826, 829, 829 
# Resampling results across tuning parameters:
# 
#   max_depth  gamma  min_child_weight  Accuracy   Kappa      LogLoss 
#   5           1     2                 0.5877763  0.4769229  1.229919
#   5           1     4                 0.5684582  0.4500727  1.249492
#   5           1     8                 0.5588254  0.4350412  1.280321
#   5           4     2                 0.5676162  0.4435441  1.311165
#   5           4     4                 0.5587963  0.4314523  1.327369
#   5           4     8                 0.5628259  0.4349619  1.360671
#   5          12     2                 0.5048484  0.3458072  1.535380
#   5          12     4                 0.5040646  0.3500563  1.544040
#   5          12     8                 0.4903788  0.3312398  1.564051
#   7           1     2                 0.5893963  0.4793720  1.236464
#   7           1     4                 0.5619897  0.4436735  1.247619
#   7           1     8                 0.5588021  0.4365085  1.287118
#   7           4     2                 0.5660311  0.4419835  1.313589
#   7           4     4                 0.5579600  0.4303822  1.337724
#   7           4     8                 0.5442568  0.4114695  1.363229
#   7          12     2                 0.5201717  0.3711661  1.537965
#   7          12     4                 0.5209905  0.3669324  1.550363
#   7          12     8                 0.4976544  0.3355580  1.572314
#   9           1     2                 0.5805589  0.4674446  1.225970 (!)
#   9           1     4                 0.5781318  0.4626184  1.244306
#   9           1     8                 0.5612292  0.4387789  1.289420
#   9           4     2                 0.5579600  0.4324878  1.315623
#   9           4     4                 0.5571646  0.4301692  1.321481
#   9           4     8                 0.5523627  0.4229930  1.362065
#   9          12     2                 0.5233827  0.3769908  1.533696
#   9          12     4                 0.5137440  0.3602794  1.537083
#   9          12     8                 0.4952098  0.3343036  1.577355
# 
# Tuning parameter 'nrounds' was held constant at a value of 75
# Tuning parameter 'eta' was held constant at a
#  value of 0.07
# Tuning parameter 'colsample_bytree' was held constant at a value of 0.1
# Tuning parameter
#  'subsample' was held constant at a value of 0.5
# LogLoss was used to select the optimal model using  the smallest value.
# The final values used for the model were nrounds = 75, max_depth = 9, eta = 0.07, gamma = 1, colsample_bytree
#  = 0.1, min_child_weight = 2 and subsample = 0.5.

# nci_gene_dtm ------------------------------------------------------------

# eXtreme Gradient Boosting 
# 
# 1242 samples
#  941 predictor
#    9 classes: 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9' 
# 
# No pre-processing
# Resampling: Cross-Validated (3 fold) 
# Summary of sample sizes: 827, 830, 827 
# Resampling results across tuning parameters:
# 
#   max_depth  gamma  min_child_weight  Accuracy   Kappa      LogLoss 
#   5           1     2                 0.5691874  0.4488836  1.263252
#   5           1     4                 0.5667602  0.4443424  1.286545
#   5           1     8                 0.5450560  0.4136148  1.325501
#   5           4     2                 0.5514758  0.4214118  1.351541
#   5           4     4                 0.5595489  0.4311033  1.362643
#   5           4     8                 0.5410750  0.4055695  1.397905
#   5          12     2                 0.5023804  0.3458818  1.571904
#   5          12     4                 0.4927243  0.3294723  1.572560
#   5          12     8                 0.4887082  0.3271582  1.585083
#   7           1     2                 0.5740476  0.4565519  1.258615
#   7           1     4                 0.5700140  0.4494941  1.284124
#   7           1     8                 0.5498869  0.4214516  1.316762
#   7           4     2                 0.5675927  0.4420510  1.343262
#   7           4     4                 0.5611904  0.4311549  1.364629
#   7           4     8                 0.5466507  0.4137830  1.392851
#   7          12     2                 0.4911179  0.3252532  1.559927
#   7          12     4                 0.5104652  0.3543722  1.561541
#   7          12     8                 0.4790755  0.3117119  1.590391
#   9           1     2                 0.5700082  0.4526958  1.263141
#   9           1     4                 0.5587691  0.4358111  1.282771
#   9           1     8                 0.5506901  0.4239392  1.328893
#   9           4     2                 0.5563126  0.4252951  1.350641
#   9           4     4                 0.5491188  0.4156286  1.356915
#   9           4     8                 0.5410340  0.4048539  1.397687
#   9          12     2                 0.5185382  0.3632035  1.557768
#   9          12     4                 0.4911530  0.3299017  1.561991
#   9          12     8                 0.4661598  0.2945271  1.588648
# 
# Tuning parameter 'nrounds' was held constant at a value of 75
# Tuning parameter 'eta' was held constant at a
#  value of 0.07
# Tuning parameter 'colsample_bytree' was held constant at a value of 0.1
# Tuning parameter
#  'subsample' was held constant at a value of 0.5
# LogLoss was used to select the optimal model using  the smallest value.
# The final values used for the model were nrounds = 75, max_depth = 7, eta = 0.07, gamma = 1, colsample_bytree
#  = 0.1, min_child_weight = 2 and subsample = 0.5.
