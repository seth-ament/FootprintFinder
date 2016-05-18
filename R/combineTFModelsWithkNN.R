# cluster motifs by footprints
# test performance of classifier using the k=5 most similar TFs

require( caret )
require( pROC )
require( doParallel )
require( foreach )
registerDoParallel( cores = 8 )


# for which TFs do we have feature matrices?
feature_matrices = Sys.glob("FeatureMatrix/*RData")
tflist = gsub("FeatureMatrix\\/" , "" , feature_matrices )
tflist = gsub("\\_(.*)" , "" , tflist )

feature_summary = matrix( NA , ncol = 39 , nrow = length(tflist) )
for( i in 1:length(tflist) ) {
  cat( i , "\n" )
  tf = tflist[i]
  featurefile = paste( "FeatureMatrix/" , tf , "_feature_matrix.RData" , sep="" )
  load( featurefile )
  x = features[,-c(1:4)]
  s1 = quantile( x[,1] , seq( 0 , 1 , 0.25 ) )
  s2 = quantile( x[,2] , seq( 0 , 1 , 0.25 ) )
  s3.30 = colSums( x[,3:30] ) / nrow(x)
  feature_summary[i,] = c( nrow(x) , s1 , s2 , s3.30 )
}

for( j in 1:ncol(feature_summary) ) {

means = colMeans( feature_summary )
sd = apply( feature_summary , 2 , sd )
norm = sapply( 1:ncol(feature_summary) , function(x) feature_summary[,x] - means

norm =  t( ( t(feature_summary) - means ) / sd )
rownames(norm) = tflist

d = as.matrix(dist( norm , upper = T ))
rownames(d) = colnames(d) = tflist
k = 10
tflist[ order( d[,3] )[1:k] ]

# load the single TF models
modelfiles = Sys.glob("SingleTFModels/*.RData")
models = list()
for( m in modelfiles ) {
  tf.m = gsub( "SingleTFModels\\/" , "" , m )
  tf.m = gsub( "\\_(.*)" , "" , tf.m )
  load( m )
  models[[tf.m]] = fit
}

# predict TFBSs in test set
# majority vote of all other models
performance.metrics = list()
foreach( tf=tflist ) %dopar% {
  cat( "working on" , tf , "\n" )
  featurefile = paste( "FeatureMatrix/" , tf , "_feature_matrix.RData" , sep="" )
  load( featurefile )
  y = features$y
  x = features[,-c(1:4)]
    kNN10.tfs = tflist[ order( d[,tf] )[1:6] ]
  predictions = matrix( NA , nrow = nrow(x) , ncol = length(kNN10.tfs) )
  colnames(predictions) = kNN10.tfs
  roc.single = rep( NA , length(kNN10.tfs) )
  for( i in 1:length(kNN10.tfs) ) {
    cat( i , "\n" )
    m = which( tflist %in% kNN10.tfs )[i]
    probs = predict( models[[m]] , newdata = x , type = "prob" )
    predictions[,i] = probs[,2]
    roc.single[i] = as.numeric( roc( response = y , predictor = probs[,2] , direction = "<" )$auc )
  }
  median.probs = apply( predictions[ , -which( kNN10.tfs == tf ) ] , 1 , quantile , probs = 0.5 )
  t2 = table( y == "yesChIP" , predictions[,tf] > 0.5 )
  sensitivity.train = t2[2,2] / sum(t2[2,])
  specificity.train = t2[2,2] / sum(t2[,2])
  t = table( y == "yesChIP" , median.probs > 0.5 )
  sensitivity.median = t[2,2] / sum(t[2,])
  specificity.median = t[2,2] / sum(t[,2])
  roc.median = roc( response = y , predictor = median.probs , direction = "<" )
  roc.training = roc( response = y , predictor = predictions[ , tf ] )
  roc.fimo = roc( response = y , predictor = x$min.p.fimo )
  roc.wellington = roc( response = y , predictor = x$max.score.wellington )
  performance.metrics = list(
        probs.singletf = predictions ,
        auc.singletf = roc.single ,
        probs.median = median.probs ,
        confusion.matrix.median = t ,
        confusion.matrix.train = t2 ,
        auc.train = as.numeric( roc.training$auc ) ,
        auc.median = as.numeric( roc.median$auc ) ,
        auc.fimo = as.numeric( roc.fimo$auc ) ,
        auc.wellington = as.numeric( roc.wellington$auc ) ,
        sensitivity.train = sensitivity.train ,
        specificity.train = specificity.train ,
        sensitivtiy.median = sensitivity.median ,
        specificity.median = specificity.median )
  outfile = paste( "PerformanceMetrics/" , tf , "_performance_metrics.RData" , sep = "" )
  save( performance.metrics , file = outfile )
}






pdf("roc.median_probs.knn5.pdf")
  plot( roc.training , lty = 2 )
  par( new = T )
  plot( roc.median , lty = 1 , col = "red" , lwd = 3 )
  par( new = T )
  plot( roc.fimo , lty = 3 )
  par( new = T )
  plot( roc.wellington , lty = 4 )
  legend( x = 0.4 , y = 0.2 ,
        lty = c(1,2,3,4) ,
        col = c("red",rep("black",3)) ,
        lwd = c(3,1,1,1) ,
        legend = c( "test (ensemble)" , "training" , "FIMO only" , "Wellington only" )
        )
  mtext( side = 3 , line = 2 , adj = 0 , font = 4 , cex = 2 , tf )
dev.off()








