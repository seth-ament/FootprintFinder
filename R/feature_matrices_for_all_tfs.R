# create feature matrices for all 850 TFs

require(TReNA)
require(doParallel)
require(foreach)
registerDoParallel( cores = 8 )


source("/proj/price1/sament/lymphoblast_trn/FootprintFinder/R/assembleFeatureMatrixNoChIP.R")

   genome.db.uri <- "postgres://whovian/hg38"
   project.db.uri <-  "postgres://whovian/lymphoblast"
   fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

query <-
paste( "select distinct tf from motifsgenes" )
 
all.tfs = dbGetQuery( fp@project.db , query )[,1]

closeDatabaseConnections(fp)

setwd("/proj/price1/sament/lymphoblast_trn")

foreach( tf=all.tfs ) %dopar% {
  cat( "###### Working on" , tf , "##############\n" )
  features = assembleFeatureMatrix( tf )
  filename = paste( "AllFeatureMatrix/" , tf , ".RData" , sep="" )
  save( features , file=filename )
}




