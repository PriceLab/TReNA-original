trimTRN <-
function( trn , r2.threshold = 0.5 , beta.threshold = 0.01 ) {

r2 = trn$r2.predVactual
coef = trn$beta.coefficients

#select genes with accurate predictions
goodPred = which( r2 >= r2.threshold )
tmp = paste( "Removed" , nrow(coef) - length(goodPred) , "genes that were not accurately predicted (R^2 <" , r2.threshold , ")\n" )
cat(tmp)
trn2 = coef[ goodPred , ]

#set edges with very small effects to zero
small.coef = which( abs(trn2) < beta.threshold & trn2 != 0 )
trn2[ small.coef ] = 0
tmp2 = paste( "Removed" , length(small.coef) , "edges with small effect sizes (beta <" , beta.threshold , ")\n" )
cat( tmp2 )

#remove TFs that do not have any predicted targets in the trimmed network
#this also removes the "intercept" column when the TRN uses zero-centered expression data.
x = which( colSums( trn2 ) == 0 )
trn2 = trn2[ , -x ]
tmp3 = paste( "Removed" , length(x) , "TFs with no predicted targets.\n" )
cat( tmp3 )

nEdges = length( which( trn2 != 0 ) )
nReg = ncol( trn2 )
nTargets = nrow( trn2 )

tmp4 = paste( "Trimmed TRN model contains" , nTargets , "target genes and" , nReg , "TFs, linked by" , nEdges , "edges.\n" )
cat( tmp4 )

trn2

}
