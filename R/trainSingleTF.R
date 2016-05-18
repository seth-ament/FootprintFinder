# train models for each TF

require( caret )
require( FootprintFinder )
require( doMC )
require( pROC )
registerDoMC( cores = 10 )

source("/proj/price1/sament/lymphoblast_trn/FootprintFinder/R/assembleFeatureMatrix.R")


tflist = Sys.glob("/proj/price1/sament/lymphoblast_trn/known_tfbs/hg38/*.bed")
tflist = gsub("/proj/price1/sament/lymphoblast_trn/known_tfbs/hg38/","",tflist)
tflist = gsub("\\_(.*)" , "" , tflist )

for( tf in tflist ) {

donetfs = Sys.glob("/proj/price1/sament/lymphoblast_trn/SingleTFModels/*.RData")
donetfs = gsub("/proj/price1/sament/lymphoblast_trn/SingleTFModels/","",donetfs)
donetfs = gsub("_gbm_model.RData","",donetfs)

if( tf %in% donetfs ) next

cat( "############Working on" , tf , "###########\n" )

cat("Assembling feature matrix\n")
features = assembleFeatureMatrix( tf )

if( is.null(features) ) next

save( features ,  file=paste( "FeatureMatrix/" , tf , "_feature_matrix.RData" , sep = "" ))

y = features$y

x = features[,-c(1:4)]

cat("fitting model\n")
ctrl = trainControl( method = "repeatedcv" ,
        repeats = 3 ,
        classProbs = TRUE )
grid = expand.grid( n.trees = c(500) ,
		interaction.depth = c(8) ,
		shrinkage = 0.1 ,
		n.minobsinnode = 10 )
fit = train( x = x , y = y ,
        method = "gbm" ,
        trControl = ctrl ,
	metric = "Kappa" ,
	tuneGrid = grid ,
	verbose = F )

save( fit , file=paste("SingleTFModels/",tf,"_gbm_model.RData" , sep="" ) )

}

#probs = predict( fit , newdata = x.test , type = "prob" )
#roc = roc( response = y.test , predictor = probs[,2] )

#table( y.test , probs[,2] > 0.5 )




