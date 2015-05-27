summarystats <-
function( trn ) {

if( is.list( trn ) ) {
	trn = trn$beta.coefficients[,-1]
}

kIn = rowSums( trn != 0 )
kIn.median = median( kIn )
kOut = colSums( trn != 0 )
kOut.median = median( kOut )
nEdges = length( which( trn != 0 ))
nTargets = nrow(trn) 
nRegulators = ncol(trn)

list( nTargets = nTargets , nRegulators = nRegulators , nEdges = nEdges , kIn.median = kIn.median , kOut.median = kOut.median )

}
