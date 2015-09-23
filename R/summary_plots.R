summary_plots =
function( trn ) {

trn.trim = trimTRN( trn )

layout( mat = matrix( c( 1 , 1 , 2 , 3 ) , nrow = 2 , ncol = 2 , byrow = T ) )

# plot of prediction accuracies
r2 = trn$r2.predVactual
na = which( is.na( r2 ) )
par( bty = "l" )
plot( 	y = sort( r2[ -na ] , decreasing = T ) , 
		x = 1:length(r2[-na]) ,
		xlab = "Genes" ,
		ylab = "Prediction Accuracy (R^2)" ,
		cex = 2 , type = "l" , lwd = 2
		)
mtext( side = 3 , "Training Set Prediction Accuracy for each Target Gene" )
n = length( which( r2 > 0.5 ) )
abline( v = n )


# plot of in-degree distribution
hist( 	rowSums( trn.trim != 0 ) , 
		xlab = "in-degree of target genes" , 
		main = "in-degree of target genes" ,
		col = "grey" , breaks = 50 , cex = 1.5
		)

# plot of out-degree distribution
hist( 	colSums( trn.trim != 0 ) , 
		xlab = "out-degree of TFs" , 
		main = "out-degree of TFs" ,
		col = "grey" , breaks = 50 , cex = 1.5 
		)


}
