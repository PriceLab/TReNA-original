options(stringsAsFactors=F)

#get all the CHDI mouse cortex read counts

###################
# HDAC DATA INPUT #
###################

#all read count files
#hdac
setwd("/proj/price1/CHDI/hdac/analysis")

load("deg/r62.crosses.gene_name.counts.RData")
counts.r62.crosses = counts
meta.r62.crosses = read.csv("/proj/price1/CHDI/hdac/metadata/Decoder_ring_R62_HDAC4_Crosses.csv")

load("deg/r62.hdaci.gene_name.counts.RData")
counts.r62.hdaci = counts
meta.r62.hdaci = read.csv("/proj/price1/CHDI/hdac/metadata/Decoder_ring_R62_HDACi_vehicle_treated.csv")

load("deg/q175.crosses.gene_name.counts.RData")
counts.q175.crosses = counts
meta.q175.crosses = read.csv("/proj/price1/CHDI/hdac/metadata/Decoder_ring_Q175_HDAC4_Crosses.csv")

load("deg/q175.hdaci.gene_name.counts.RData")
counts.q175.hdaci = counts
meta.q175.hdaci = read.csv("/proj/price1/CHDI/hdac/metadata/Decoder_ring_Q175_HDACi_vehicle_treated.csv")

rm(counts)
tmp1 = which( meta.r62.crosses$Tissue == "CORTEX" ) #99
tmp2 = grep( "cortex" , meta.r62.hdaci$Source2 ) #32
tmp3 = which( meta.q175.crosses$Tissue == "CORTEX" ) #18
tmp4 = grep( "cortex" , meta.q175.hdaci$Tissue ) #32

counts = cbind( counts.r62.crosses[ , tmp1 ] , counts.r62.hdaci[ , tmp2 ] , counts.q175.crosses[ , tmp3 ] , counts.q175.hdaci[ , tmp4 ] )

tmp1 = which( meta.r62.crosses$Tissue == "CORTEX" ) #99
tmp2 = grep( "cortex" , meta.r62.hdaci$Source2 ) #32
tmp3 = which( meta.q175.crosses$Tissue == "CORTEX" ) #18
tmp4 = grep( "cortex" , meta.q175.hdaci$Tissue ) #32

counts.hdac = cbind( counts.r62.crosses[ , tmp1 ] , counts.r62.hdaci[ , tmp2 ] , counts.q175.crosses[ , tmp3 ] , counts.q175.hdaci[ , tmp4 ] )

meta1 = cbind( meta.q175.crosses , NA  , "hdac.q175.crosses" )
meta2 = cbind(  meta.q175.hdaci[,1:2] , NA , meta.q175.hdaci[,c(4,3,5:6)] , "hdac.q175.hdaci" )
meta3 = cbind(  meta.r62.crosses , NA  , "hdac.r62.crosses" )
meta4 = cbind(  meta.r62.hdaci[,1] , meta.r62.hdaci[,7] ,NA , meta.r62.hdaci[,c(4,3,5,6)] , "hdac.r62.hdaci" )
names(meta3) = names(meta4) = names(meta1) = names(meta2) = c(  names(meta1)[1:6] , "Sex" , "Study" )
meta.hdac = rbind( meta1 , meta2 , meta3 , meta4 )

meta.hdac$Tissue = tolower( meta.hdac$Tissue )
meta.hdac$Tissue = gsub( "l\\." , "" , meta.hdac$Tissue )
meta.hdac$Tissue = gsub( "left_" , "" , meta.hdac$Tissue )
meta.hdac$Tissue = gsub( "tibilis" , "tibialis" , meta.hdac$Tissue )

matchLines = match( colnames(counts.hdac) , meta.hdac[,1] )
meta.hdac = meta.hdac[ matchLines , ]


#############################
# ALLELIC SERIES DATA INPUT #
#############################

setwd("/proj/price1/CHDI/Allelic_series/metadata2")
expt = "allelic_series"
metadata = read.csv("allelic_series_master_decoder_ring.csv")
useSamples = which(metadata$Sequence %in% c("mRNA") )
mRNAdata = metadata[useSamples,]
meta.samples = substr( mRNAdata , 1 , 6 )

#locations of all the count files from individual samples
setwd("/proj/price1/jpearl/snapr/Allelic_series")
proc = Sys.glob("./*/*/*gene_name.counts.txt")
x = strsplit( proc , split = "/" )
proc.id = rep( NA , length( x ) ) 
for( i in 1:length(x) ) {
	proc.id[i] = x[[i]][3]
}
trim.id = substr( proc.id , 1 , 6 )

#create a count table for all the samples in this experiment.
x = read.table( proc[1] )
nSamples = length(proc)
nGenes = nrow(x)
counts = matrix( NA , nGenes , nSamples )
rownames(counts) = x[,1]
colnames(counts) = trim.id
for( i in 1:nSamples ) {
	sample.counts = read.table( proc[i] )
	counts[ , i ] = sample.counts[,2]
	cat(i)
	cat("\n")
}

noData = which( colSums( counts ) == 0 )
counts = counts[ , -noData ]

meta = cbind( meta.samples , mRNAdata )
rownames(meta) = 1:nrow(meta)
matchLines = match( colnames(counts) , meta.samples )
meta = meta[ matchLines , ]
meta = meta[ , 1:14 ]

meta5 = cbind( meta[ , c(1,8,11) ] , NA , meta[,c(5,4,7)] , "allelic_series" )
names( meta5 ) = names( meta.hdac )

cortex = which( meta$Tissue == "cortex" )
counts.allelic_series = counts[,cortex]

counts = cbind( counts.hdac , counts.allelic_series )

meta.combined = rbind( meta.hdac , meta5 )
matchLines = match( colnames(counts) , meta.combined[,1] )
meta.combined = meta.combined[ matchLines , ]

setwd("/proj/price1/CHDI/users/sament")
save( counts , file = "combined_cortex_read_counts.RData" )
save( meta.combined , file = "combined_cortex_metadata.RData" )

conds = meta.combined

library(DESeq)
cds <- newCountDataSet( counts , conds )
cds <- estimateSizeFactors( cds )
norm.counts = t( t(counts(cds)) / sizeFactors(cds) )

library( WGCNA )









