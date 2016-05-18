# Cluster the TFs based on features of their footprints

# for which TFs do we have feature matrices?
feature_matrices = Sys.glob("AllFeatureMatrix/*RData")
tflist = gsub("AllFeatureMatrix\\/" , "" , feature_matrices )
tflist = gsub(".RData" , "" , tflist )

feature_summary = matrix( NA , ncol = 39 , nrow = length(tflist) )
for( i in 1:length(tflist) ) {
  cat( i , "\n" )
  tf = tflist[i]
  featurefile = paste( "AllFeatureMatrix/" , tf , ".RData" , sep="" )
  load( featurefile )
  x = features[,-c(1:3)]
  s1 = quantile( x$min.p.fimo , seq( 0 , 1 , 0.25 ) )
  s2 = quantile( x$max.score.wellington , seq( 0 , 1 , 0.25 ) )
  s3.30 = colSums( x[,-c(1:2)] ) / nrow(x)
  feature_summary[i,] = c( nrow(x) , s1 , s2 , s3.30 )
}



means = colMeans( feature_summary )
sd = apply( feature_summary , 2 , sd )
norm =  t( ( t(feature_summary) - means ) / sd )
rownames(norm) = tflist

d = as.matrix(dist( norm , upper = T ))
rownames(d) = colnames(d) = tflist

save( d , file="footprint_features_distance_matrix.RData" )




