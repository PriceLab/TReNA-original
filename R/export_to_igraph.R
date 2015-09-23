trena2igraph =
function( trn , cols = NULL, rows = NULL , beta.threshold = 0 ) {
library(igraph)
if( is.null(cols) ) cols = colnames(trn)
if( is.null(rows) ) rows = rownames(trn)
trn = trn[ which(rownames(trn) %in% rows) , which(colnames(trn) %in% cols) ]
x = which( abs(trn) > beta.threshold , arr.ind = T )
source = colnames(trn)[ x[,2] ]
target = rownames(trn)[ x[,1] ]
weight = abs (trn[x] )
beta = trn[x]
sign = sign(trn[x])
graph = data.frame( source , target , beta , weight , sign )
graph_from_data_frame( graph , directed = T )
}
