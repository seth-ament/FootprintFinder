assembleFeatureMatrix = function( tf ) {

require( GenomicRanges )
require( TReNA )
require( doBy )

cat( "initializing databases\n" )
genome.db.uri <- "postgres://whovian/hg38"
project.db.uri <-  "postgres://whovian/lymphoblast"
fp <- FootprintFinder(genome.db.uri, project.db.uri, quiet=TRUE)

###### Wellington / FIMO ##########

cat( "loading footprints\n" )
# get Wellington/FIMO footprints for tf

query <-
paste( "select fp.chr, fp.mfpstart, fp.mfpend, fp.motifname, fp.pval, fp.score, mg.motif, mg.tf",
             "from footprints fp",
             "inner join motifsgenes mg",
             "on fp.motifName=mg.motif",
             paste("where mg.tf='",tf,"'" , sep="" ) ,  collapse = " " )

footprints.tf = dbGetQuery( fp@project.db , query )

colnames( footprints.tf ) = c("chrom","start","end","motifname","p.fimo","score.wellington","motif","tf" )
fp.gr = makeGRangesFromDataFrame( footprints.tf , keep.extra.columns = T )
footprints.uniq = unique(footprints.tf[,1:3])
fp.gr.uniq = makeGRangesFromDataFrame( footprints.uniq )

cat( "storing FIMO and Wellington p-values" )

# min FIMO p-value

loc = paste( footprints.tf$chr , footprints.tf$start ,
        footprints.tf$end , footprints.tf$tf , sep="." )

x = data.frame( loc , footprints.tf )
x = x[ order( x$p.fimo ) , ]
x = x[ duplicated(x$loc) == F , ]
m = match( unique(loc) , x$loc )
min.p.fimo = x$p.fimo[ m ]

# max Wellington score

x = data.frame( loc , footprints.tf )
x = x[ order( x$score.wellington , decreasing = T ) , ]
x = x[ duplicated(x$loc) == F , ]
m = match( unique(loc) , x$loc )
max.score.wellington = x$score.wellington[ m ]


###### ChIP-seq (gold standard) #######

cat( "loading ChIP-seq\n" )
# ChIP-seq peaks for tf (in lymphoblasts)
chipfile = paste( "/proj/price1/sament/lymphoblast_trn/known_tfbs/hg38_train/", tf , "_lymphoblast_binding_sites.bed.hg38.bed" , sep="" )
chip.tf = read.table(chipfile)[,1:4]
colnames(chip.tf) = c("chrom","start","end","tfs")
chip.gr = makeGRangesFromDataFrame( chip.tf )
# matches to ChIP (true/false positives)
hasChIP = countOverlaps( fp.gr.uniq , chip.gr , maxgap = 100 )
hasChIP = hasChIP > 0

# y is the outcome variable for prediction
y = rep(NA,length(hasChIP))
y[ hasChIP == TRUE ] = "yesChIP"
y[ hasChIP == FALSE ] = "noChIP"
y = factor(y)

####### phastCons ########

cat( "loading phastCons\n" )
# get phastCons conserved element annotations
load("/proj/price1/sament/resources/phastConsElements100way.RData")
phastCons.gr = makeGRangesFromDataFrame( phastCons )
# matches to conserved elements
conserved_element = countOverlaps( fp.gr.uniq , phastCons.gr )

###### ChromHMM ########

# imputed 25 state model
# get chromHMM annotations for GM12878

cat( "loading ChromHMM\n" )
hmm.states = read.csv("/proj/price1/sament/resources/ChromHMM/states.25statemodel.csv" , header = F)
chromHMM = read.table( "/proj/price1/sament/lymphoblast_trn/E116_25_imputed12marks_stateno.hg38.bed")
colnames(chromHMM) = c("chr","start","end","state")
chromHMM = chromHMM[ order( chromHMM[,4] ) , ]
chromStates = splitBy( ~ state , data = chromHMM )

# matches to each chromHMM state

fp.chromStates = matrix( NA , 
	ncol = length(chromStates) , 
	nrow = nrow(footprints.uniq) )

colnames(fp.chromStates) = paste( "chromHMM" , hmm.states[,2] , sep="." )

for( i in 1:length(chromStates) ) {
  df = chromStates[[i]]
  chromState.gr = makeGRangesFromDataFrame( df )
  counts.state = countOverlaps( fp.gr.uniq , chromState.gr )
  fp.chromStates[,i] = counts.state
}

####### FANTOM5 #########

cat( "loading FANTOM5\n" )

# FANTOM5 nehancers B-lymphocytes
fantom = read.table( "/proj/price1/sament/resources/fantom5/hg38/facet_expressed_enhancers/CL0000945_lymphocyte_of_B_lineage_expressed_enhancers.bed.hg38.bed")
colnames(fantom)[1:3] = c("chrom","start","end")
fantom = makeGRangesFromDataFrame( fantom )
fantom = countOverlaps( fp.gr.uniq , fantom )

# FANTOM5 enhancers, all tissues
fantom.files = Sys.glob("/proj/price1/sament/resources/fantom5/hg38/facet_expressed_enhancers/*.bed")

fantom.all = read.table(fantom.files[1])[,1:3]
for( f in fantom.files[-1] ) {
  tmp = read.table(f)[,1:3]
  fantom.all = rbind( fantom.all , tmp )
  # cat( f , "\n" )
}
colnames(fantom.all) = c("chrom","start","end")
fantom.all = makeGRangesFromDataFrame( fantom.all )
fantom.all = countOverlaps( fp.gr.uniq , fantom.all )

# assemble feature matrix
cat( "assembling feature matrix\n" )
features = data.frame( 
	footprints.uniq ,
	y ,
	min.p.fimo , 
	max.score.wellington ,
	fp.chromStates , 
	fantom.lymphocytes = fantom ,
	fantom.all ,
	conserved_element )

return( features )

}

