# combine single-TF FootprintFinder models to predict TF occupancy for other TFs
# Seth Ament
# 2016-5-17

require( caret )
require( pROC )
require( gbm )

# for which TFs do we have feature matrices?
feature_matrices = Sys.glob("FeatureMatrix/*RData")
tflist = gsub("FeatureMatrix\\/" , "" , feature_matrices )
tflist = gsub("\\_(.*)" , "" , tflist )

# load the single TF models
modelfiles = Sys.glob("SingleTFModels/*.RData")
models = list()
for( m in modelfiles ) {
  tf.m = gsub( "SingleTFModels\\/" , "" , m )
  tf.m = gsub( "\\_(.*)" , "" , tf.m )
  load( m )
  models[[tf.m]] = fit
}

# compare the variable importance scores across models
nmodels = length(models)
nvars = nrow( as.data.frame( summary( models[[1]] ) ) )
varImp = matrix( NA , ncol = nmodels , nrow = nvars )
for( m in 1:length(models) ) {
  rel.inf = as.data.frame( summary( models[[m]] ) )
  rel.inf = rel.inf[ order( rel.inf[,1] ) , ]
  rel.inf = rel.inf[,2] / max(rel.inf[,2])
  varImp[ , m ] = rel.inf
}
colnames(varImp) = names(models)
rownames(varImp) = sort( as.data.frame( summary( models[[m]] ) )[,1] )
# heatmap to visualize
my_palette = colorRampPalette(c("white","blue"))(299)
pdf("varImp.heatmap.43singleTFModels.pdf")
heatmap(varImp , 
	col = my_palette ,
	scale = "none" ,
	margins = c(10,15) )
dev.off()
write.csv( varImp , file="variable.importance.43tfs.csv" )

# predict TFBSs in test set
# majority vote of all other models
pdf("roc.median_probs_43tfs.pdf")
performance.metrics = list()
for( tf in tflist ) {
  cat( "working on" , tf , "\n" )
  featurefile = paste( "FeatureMatrix/" , tf , "_feature_matrix.RData" , sep="" )
  load( featurefile )
  y = features$y
  x = features[,-c(1:4)]
  predictions = matrix( NA , nrow = nrow(x) , ncol = nmodels )
  colnames(predictions) = names(models)
  roc.single = rep( NA , nmodels )
  for( m in 1:nmodels ) { 
    probs = predict( models[[m]] , newdata = x , type = "prob" )
    predictions[,m] = probs[,2]
    roc.single[m] = as.numeric( roc( response = y , predictor = probs[,2] , direction = "<" )$auc )
  }
  median.probs = apply( predictions[ , -which( names(models) == tf ) ] , 1 , quantile , probs = 0.5 )
  t2 = table( y == "yesChIP" , predictions[,names(models) == tf] > 0.5 )
  sensitivity.train = t2[2,2] / sum(t2[2,])
  specificity.train = t2[2,2] / sum(t2[,2])
  t = table( y == "yesChIP" , median.probs > 0.5 )
  sensitivity.median = t[2,2] / sum(t[2,])
  specificity.median = t[2,2] / sum(t[,2])
  roc.median = roc( response = y , predictor = median.probs , direction = "<" )
  roc.training = roc( response = y , predictor = predictions[ , names(models) == tf ] )
  roc.fimo = roc( response = y , predictor = x$min.p.fimo )
  roc.wellington = roc( response = y , predictor = x$max.score.wellington )
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
  performance.metrics[[tf]] = list( 
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
}
dev.off()


quintiles.batf = t( sapply( 1:nrow(predictions) ,
        function(x) quantile( predictions[ x ,-which( names(models)==tf)] , probs = seq(0,1,0.2) )))


require(doMC)
registerDoMC( cores = 5 )
ctrl = trainControl( method = "repeatedcv" ,
		repeats = 2 )
fit.combined = train( y = y , x = quintiles , 
	trControl = ctrl , method = "gbm" , tuneLength = 2 , metric = "Kappa" )
prob.combined = predict( fit.combined , newdata = quintiles , type = "prob" )
plot( roc( response = y , predictor = prob.combined[,2] ) , lty = 2 , col = "green" )

table( y , prob.combined[,2] > 0.5 )

pred = prediction( prob.combined[,2] , y )
perf = performance( pred , measure="prec", x.measure="rec" )
pdf("ikzf1.precrecall.pdf")
plot( performance( prediction( 
	predictions[,names(models) == tf] , y ) , 
	measure="prec", x.measure="rec" ) , lty = 1 , ylim = c(0,1) )
par( new = T )
plot( performance( prediction( prob.combined[,2] , y ) , measure="prec", x.measure="rec" ) , 
	lty =2 , ylim = c(0,1) )
par( new = T )
plot( performance( prediction( median.probs , y ) , measure="prec", x.measure="rec" ) , 
        lty =2 , ylim = c(0,1) )
dev.off()








