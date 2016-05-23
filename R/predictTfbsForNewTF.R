# cluster motifs by footprints
# test performance of classifier using the k=5 most similar TFs

require( caret )
require( doParallel )
require( foreach )
registerDoParallel( cores = 4 )


# for which TFs do we have feature matrices?
feature_matrices = Sys.glob("AllFeatureMatrix/*RData")
tflist = gsub("AllFeatureMatrix\\/" , "" , feature_matrices )
tflist = gsub(".RData" , "" , tflist )

# load the single TF models for ensemble prediction
modelfiles = Sys.glob("SingleTFModels/*.RData")
models = list()
for( m in modelfiles ) {
  tf.m = gsub( "SingleTFModels\\/" , "" , m )
  tf.m = gsub( "\\_(.*)" , "" , tf.m )
  load( m )
  models[[tf.m]] = fit
}

# load the TFBS distance matrix
load("footprint_features_distance_matrix.RData")
d = d[ rownames(d) %in% names(models) , ]

donetfs = Sys.glob("FittedTFBSProbabilities/*")
donetfs = gsub("FittedTFBSProbabilities/","",donetfs)
donetfs = gsub(".RData","",donetfs)

tflist = setdiff( tflist , donetfs )

# predict TFBSs in test set
# majority vote of k models
k = 10
foreach( tf=tflist ) %dopar% {
  cat( "working on" , tf , "\n" )
  featurefile = paste( "AllFeatureMatrix/" , tf , ".RData" , sep="" )
  load( featurefile )
  x = features[,-c(1:3)]
  if( tf %in% rownames(d) ) {
    d.sub = d[ -which( rownames(d) == tf) , ]
  } else d.sub = d
  kNN = rownames(d.sub)[ order( d.sub[,tf] )[1:k] ]
  predictions = matrix( NA , nrow = nrow(x) , ncol = k )
  colnames(predictions) = kNN
  for( i in 1:k ) {
    #cat( i , "\n" )
    m = kNN[i]
    probs = predict( models[[m]] , newdata = x , type = "prob" )
    predictions[,i] = probs[,2]
  }
  median.probs = apply( predictions , 1 , median )
  median.scaled = median.probs / max(median.probs)
  mean.prob = rowMeans( predictions )
  mean.scaled = mean.prob / max(mean.prob)
  bed = cbind( features[,1:3] , median.probs , median.scaled , mean.prob , mean.scaled ,predictions )
  save( bed , file = paste( "FittedTFBSProbabilities/" , tf , ".RData" , sep="" ) )
}



