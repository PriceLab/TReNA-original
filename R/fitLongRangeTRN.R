fitLongRangeTRN = 
function( expr , tfbs.proximal , tfbs.distal , lambda = NULL , alpha = 1 ) {

# Fit a TRN incorporating information about distal TFBSs

# fit a a model using only TFs with proximal binding sites
cat("fitting regression model for TFs with proximal binding sites\n")
cat("evaluating penalty parameters\n")
tmp = estLambda( expr = expr , physnet = tfbs.proximal )
trn.proximal = fitTRN( expr = expr , physnet = tfbs.proximal , lambda = tmp$lambda.median )

# update the TRN with distal binding sites
cat("fitting regression model for TFs with distal binding sites\n")
cat("evaluating penalty parameters\n")
tfbs.p2 = tfbs.proximal[ match( rownames(tfbs.distal) , rownames(tfbs.proximal) ) , ]
tfbs.p2[ is.na(tfbs.p2) ] = FALSE
tfbs.distal[ tfbs.p2 == T ] = FALSE

offset = trn.proximal$predicted.expression
offset[ is.na( offset ) ] = 0
tmp = estLambda( expr = expr , physnet = tfbs.distal , offset = offset )
trn.distal = fitTRN( expr = expr , physnet = tfbs.distal , lambda = 0.1 , offset = offset )

cat("concatenating results\n")
newbeta = trn.proximal$beta.coefficients + trn.distal$beta.coefficients

location = matrix( NA , nrow( newbeta) , ncol(newbeta) )
location[ trn.distal$beta.coefficients != 0 ] = "distal"
location[ trn.proximal$beta.coefficients != 0 ] = "proximal"

out = list( 	fits.proximal = trn.proximal$fits ,
		fits.distal = trn.distal$fits ,
		actual.expression = trn.distal$actual.expression ,
		predicted.expression = trn.distal$predicted.expression ,
		predicted.expression.proximal = trn.proximal$predicted.expression ,
		beta.coefficients = newbeta ,
		r2.predVactual = trn.distal$r2.predVactual ,
		r2.predVactual.proximal = trn.proximal$r2.predVactual ,
		location = location 
	)

return( out )

}

