options(stringsAsFactors=F)
setwd("/proj/price1/CHDI/users/sament/mouse_cortex_trn")

#load read counts and metadata
load("combined_cortex_metadata.RData")
load("combined_cortex_read_counts.RData")

#normalization with DESeq
library(DESeq)
cds <- newCountDataSet( counts , meta.combined )
cds <- estimateSizeFactors( cds )
#cds <- estimateDispersions( cds )
norm.counts = t( t(counts(cds)) / sizeFactors(cds) )

#qc with WGCNA tools
library(WGCNA)
library(flashClust)

datExpr0 = t(norm.counts)
datExpr0[ datExpr0 == 0 ] = NA
gsg = goodSamplesGenes(datExpr0, verbose = 3);
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
datExpr0[is.na(datExpr0)] = 0

cds = newCountDataSet( counts[  gsg$goodGenes , gsg$goodSamples ]  , meta.combined[ gsg$goodSamples , ] )
cds <- estimateSizeFactors( cds )

attach( meta.combined )
traitData = data.frame( as.numeric(factor(Age)) , as.numeric(factor(Genotype)) , as.numeric(factor(Background)) , as.numeric(factor(Study)) , as.numeric(factor(Sex)) , sizeFactors(cds) )
detach(meta.combined)
names(traitData) = c( "Age" , "Condition" , "Background" , "Study" , "Sex" , "sizeFactors" )
traitColors = numbers2colors( traitData )

sampleTree2 = flashClust(dist(datExpr0), method = "average")
pdf("combined_cortex_sample_heatmap.pdf")
plotDendroAndColors(sampleTree2, traitColors,
groupLabels = names(traitData),
main = "Combined Mouse Cortex Samples")
abline( h = 750000 )
dev.off()

#look at distribution of sizeFactors
q = quantile( sizeFactors(cds) )
sqrt.size = sqrt( sizeFactors(cds) )
sd = sd( sqrt.size )
mean = mean( sqrt.size )
z.score = ( sqrt.size - mean ) / sd
#there are

clust = cutreeStatic(sampleTree2, cutHeight = 750000, minSize = 10)
table(clust)

pdf( "sizeFactor_boxplots.pdf" )
boxplot( sizeFactors(cds) ~ clust )
dev.off()
#seems like clusters 0 and 2 have low read counts, as well as one outlier sample with extremely high read counts.

keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datTraits = cbind( meta.combined[ keepSamples , ] , sizeFactors(cds)[keepSamples] )

traitColors2 = traitColors[ keepSamples , ]

sampleTree = flashClust(dist(datExpr), method = "average")
pdf("combined_cortex_sample_heatmap_rm_outliers.pdf")
plotDendroAndColors(sampleTree, traitColors2,
groupLabels = names(traitData),
main = "Combined Mouse Cortex Samples (outliers removed)")
dev.off()

save( datExpr , datTraits , file = "cleaned_cortex_normalized_read_counts_and_metadata.RData")

