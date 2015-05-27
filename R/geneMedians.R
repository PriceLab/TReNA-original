geneMedians = 
function( expr , anno ) {

genes = unique( anno )
genes = setdiff( genes , c("" , NA) )
geneExpr = matrix( nrow = length(genes) , ncol = ncol(expr) )
for( i in 1:length(genes) ) {
	x = which( anno == genes[i] )
	if( length(x) == 1 ) geneExpr[i,] = expr[x,]
	if( length(x) > 1 ) {
		tmp = expr[x,]
		geneExpr[i,] = apply( tmp , 2 , median )
	}
	cat( i ) 
	cat( "\n" )
}
colnames(geneExpr) = colnames(expr)
rownames(geneExpr) = genes

geneExpr

}
