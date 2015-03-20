#Network Reconstruction
options(stringsAsFactors=F)

#load normalized, qc'ed read counts and metadata
setwd("/proj/price1/CHDI/users/sament/mouse_cortex_trn")
load("cleaned_cortex_normalized_read_counts_and_metadata.RData")

#load TFBS database
setwd("/proj/price1/sament/hdlux/")
load("tfbs_network_mmu_2014-10-23.RData")
tfnet.rows = row.names(physnet)
tfnet.rows = substr( tfnet.rows , 5 , 25 )

setwd("/proj/price1/sament")
geneinfo = read.delim("Mus_musculus.gene_info" , comment.char = "#" )

mgi = merge( colnames(datExpr) , geneinfo[ , c(2,3,5,9) ] , by.x = 1 , by.y = 2 , all.x = T , sort = T )
anno = mgi
anno = anno[ -which( duplicated( anno[,1] ) == T ) , ]
rownames( anno ) = anno[,1]

expr = t( datExpr )
expr = expr[ order( rownames( expr ) ) , ]
anno = anno[ order( anno[,1] ) , ]
has.anno = intersect( rownames( expr ) , anno[,1] )
expr = expr[ has.anno , ]
table( rownames( expr) == rownames( anno ) )

nas = which( is.na( anno$GeneID ))

expr = expr[ -nas , ]
anno = anno[ -nas , ]
table( rownames( expr) == rownames( anno ) )

actual = expr
pred = matrix( NA , nrow = nrow( actual ) , ncol = ncol( actual ) )

tf.symbols = intersect( anno[,1] , colnames( physnet ) )

library(glmnet)

fits = list(  )

r.pVa= vector(length = nrow(actual) )
r.pVa[ 1:nrow(actual) ] = NA
p.pVa = vector(length = nrow(actual) )
p.pVa[ 1:nrow(actual) ] = NA

lambda = vector(length = nrow(actual) )
lambda[ 1:nrow(actual) ] = NA
coefficients = matrix( 0 , nrow = nrow(actual) , ncol = length(tf.symbols) + 1 )
colnames(coefficients) = c("Intercept" , tf.symbols )
rownames(coefficients) = rownames(actual)
setwd("/proj/price1/CHDI/users/sament/mouse_cortex_trn")
for( i in 1:nrow(actual) ) { 
	uid = anno[i,2]
	symbol = anno[i,1]
	if( is.na(uid) ) next
	status = paste( "Working on row" , i , "," , symbol , ": " )
	cat(status)
	y = actual[ i , ]
	sub = which( tfnet.rows == uid )
	if( length(sub) == 0 ) {
		cat("\n")
		next
	}
	if( length(sub) == 1 ) {
		has_tfbs = which( physnet[sub,] > 0 )
		potential.regulators = names(has_tfbs)
	}
	if( length(sub) > 1 ) {
		has_tfbs = which( colSums(physnet[sub,]) > 0 )
		potential.regulators = names(has_tfbs)
	}
	potential.regulators = intersect( potential.regulators , tf.symbols )
	if( length(potential.regulators) == 0 ) cat("\n")
	if( length(potential.regulators) == 0 ) next
	pr = setdiff( which( anno[,1] %in% potential.regulators) , which( anno[,1]==symbol) )
	if( length(pr) < 2 ) cat("\n")
	if( length(pr) < 2 ) next
	x = t( actual[ pr , ] )
	fit = cv.glmnet( y = y , x = x , family = "gaussian" )
	fits[[i]] = fit
	s = fit$lambda.1se
	pred.i = predict( fit , newx = x , s = s , type = "link" )
	actual.i = actual[ i , ]
	r.pVa[i] = cor.test( pred.i[,1] , actual.i )$estimate
	p.pVa[i] = cor.test( pred.i[,1] , actual.i )$p.value
	pred[i,] = pred.i
	nonzero = predict( fit , newx = x , s = s , type = "nonzero" )
	lambda[i] = fit$lambda.1se
	if(!is.null(nonzero[[1]]))	cat( paste( length(nonzero[,1]) , "regulators, r =" , r.pVa[i] ) )
	coef.i = predict(  fit , newx = x , s = s , type = "coefficients" )[,1]
	names(coef.i)[1] = "Intercept"
	match.coef = match( names(coef.i) , colnames(coefficients) )
	coefficients[ i , match.coef ] = coef.i	
	cat("\n")
}
save( fits , file = "fits.RData")
save( coefficients , lambda , file = "coefficients.RData" )
save( pred , actual , p.pVa , r.pVa , file = "predictions.RData" )
save( anno , file = "rows_metadata.RData" )
