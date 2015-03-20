options(stringsAsFactors=F)

load("predictions.RData")
load("coefficients.RData")
load("rows_metadata.RData")

r2 = r.pVa ^ 2
pdf( "r2_histogram.pdf" )
par( lwd = 2 )
hist( r2 , breaks = 50 , col = "black" , cex.axis = 2 , xlab = "" , ylab = "" , main = "" , ylim = c(0,1200) )
dev.off()

goodPred = which( r2 > 0.5 )

trn = coefficients[ goodPred , ]
anno = anno[ goodPred , ]

kIn = summary ( rowSums( trn[,-1] != 0 ) ) # median = 416 )
names(kIn) = paste( "kIn" , names(kIn) , sep = "." )
kOut = summary ( colSums( trn[,-1] != 0 ) ) # median = 28 )
names(kOut) = paste( "kOut" , names(kOut) , sep = "." )
nEdges = length( which( trn[,-1] != 0 )) #393696
names(nEdges) = "nEdges"
nTargets = nrow(trn) 
names(nTargets) = "nTargets"
nRegulators = ncol(trn)-1
names(nRegulators) = "nRegulators"
stats = c( kIn , kOut , nEdges , nTargets , nRegulators )
stats = data.frame( names(stats) , stats )

write.table( stats , row.names=F , col.names=F , quote=F , sep = "\t" , file = "trn_stats.txt")

save( trn , anno , file = "trn_and_metadata.RData" )