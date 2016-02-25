fitTRN_old <-
function( expr , physnet , lambda = 0.1 , cores = NULL , alpha = 1 ) {

require(glmnet)
require(foreach)
require(doParallel)

registerDoParallel( cores = cores );
cat("doParallel initialized with ", getDoParWorkers(), " cores ", date(), "\n")

#remove columns of the TFBS (TFs) for which we do not have expression data.
physnet = physnet[ , which( colnames(physnet) %in% rownames(expr) ) ]

#center and scale the expression data.
sd = apply( expr , 1 , sd )
actual = ( expr - rowMeans(expr) ) / sd

# initialize some matrices and vectors to store information about the lasso 
# regression model - most of this is now done after foreach loop.
tf.symbols = intersect( rownames(expr) , colnames( physnet ) )
coefficient_names = matrix( 0 , nrow = 1 , ncol = length(tf.symbols) + 1 )
colnames(coefficient_names) = c("Intercept" , tf.symbols )

if( lambda != "adaptive" && !is.numeric(lambda) ) {
  cat( "lambda is out of range. please specify a numeric value or 'adaptive'." )
  return()
}

#run the lasso regression for all genes in parallel
cat("Building lasso regression models for all genes" , date() , "\n" )
loopvar<-foreach( i=1:nrow(actual) ) %dopar% {

#symbol is the ID of the target gene.
symbol = rownames(actual)[i]

#print the name of the target gene on the console
status = paste( "Working on row" , i , "," , symbol , ": " )
#cat(status)

gene_fail <- FALSE;
return_list <- list();

#y is the expression data for the target gene
y = actual[ i , ]

#sub defines the row for this target gene in the TFBS database
sub = which( rownames(physnet) == symbol )

#if the target gene is not present in the TFBS database, we skip this gene.
if( length(sub) == 0 ) {
#cat("\n")
#next
  gene_fail = TRUE
  return( return_list )
}

#most genes are represented by a single row
#has_tfbs defines the columns (TFs) that have a binding site in the promoter of this target gene.
if( length(sub) == 1 ) {
  has_tfbs = which( physnet[sub,] > 0 )
  potential.regulators = colnames(physnet)[has_tfbs]
}

#a few edge cases may have multiple rows in the TFBS database. We take the union of TFs with binding sites.
if( length(sub) > 1 ) {
  has_tfbs = which( colSums(physnet[sub,]) > 0 )
  potential.regulators = colnames(physnet)[has_tfbs]
}

#potential.regulators is a vector with the names of TFs that have binding sites in the promoter of target gene i#
#we consider only those TFs for which we have expression data, defined by the vector tf.symbols.
potential.regulators = intersect( potential.regulators , tf.symbols )

#if there are no predicted regulators, we skip this gene.
if( length(potential.regulators) == 0 ) {
  gene_fail = TRUE;
  return( return_list )
}

#autoregulation -- in which a TF regulates itself -- throws off the model. Although autoregulation is potentially interesting, we eliminate autoregulatory loops.
pr = setdiff( which( rownames(actual) %in% potential.regulators) , which( rownames(actual)==symbol) )

#if there are no predicted regulators after eliminating autoregulatory loops, we discard this gene.
if( length(pr) < 2 ) {
  gene_fail = TRUE;
  return( return_list )
}

#x is a matrix with the expression levels of TFs (potential regulators) in each of the samples.
x = t( actual[ pr , ] )

#glmnet fits a lasso regression model to predict y from the expression levels of TFs.
#there are two 'flavors', depending on how we wish to specify the penalty paramter lambda.
#if we wish to specify a uniform penalty across all target genes, we set lambda to a numeric value (default - 0.1)
#alternatively, we can run an internal cross-validation with cv.glmnet to specify an appropriate lambda for each gene.
#in the latter case, we choose the maximum value of lambda for each target gene, such that the residual variance is <1 standard deviation above the minimum residual variance (which may be overfit)

if( lambda == "adaptive" ) {
  fit = cv.glmnet( y = y , x = x , alpha = alpha , family = "gaussian" )
  s = fit$lambda.1se
} else
if( is.numeric(lambda) ) {
  fit = glmnet( y = y , x = x , alpha = alpha , family = "gaussian" )
  s = lambda
} else {
  cat( "Inconceivable!\n" )
}

#fits[[i]] = fit
return_list[[1]] = fit;

#lambda.i[i] = s
return_list[[2]] = s 

#here we predict the values of y using the regression model produced by glmnet, using the lambda specified above.
#we compare the predicted values to the actual expression levels of gene i.
#and we save the predicted values to the pred matrix.
pred.i = predict( fit , newx = x , s = s , type = "link" )
#pred[i,] = pred.i
return_list[["3"]] = pred.i

r.pVa.i = cor( pred.i[,1] , y )
return_list[["4"]] = r.pVa.i ^ 2 

#now we extract the coefficients of the TFs with non-zero estimates in our model and save these to the coefficients matrix
nonzero = predict( fit , newx = x , s = s , type = "nonzero" )
#if(!is.null(nonzero[[1]]))cat( paste( length(nonzero[,1]) , "regulators, r =" , r.pVa[i] ) )
coef.i = predict(  fit , newx = x , s = s , type = "coefficients" )[,1]
names(coef.i)[1] = "Intercept"
match.coef = match( names(coef.i) , colnames(coefficient_names) )
#coefficients[ i , match.coef ] = coef.i
return_list[[5]] = coef.i
return_list[[6]] = match.coef
#cat("\n")

return(return_list)

} # End foreach construct (do/dopar)

cat("Finished building lasso regression models for all genes, compiling results ", date(), "\n" )

#initialize some matrices and vectors to store information about the lasso 
# regression model
pred = matrix( NA , nrow = nrow( actual ) , ncol = ncol( actual ) )
fits = list()
r.pVa= rep( NA , nrow(actual) )
coefficients = matrix( 0 , nrow = nrow(actual) , ncol = length(tf.symbols) + 1 )
colnames(coefficients) = c("Intercept" , tf.symbols )
rownames(coefficients) = rownames(actual)
lambda.i = rep( NA , nrow(actual) )

# Iterate over foreach 'list of lists' to build proper TRN structure
for( i in 1:length(loopvar)) {
#  if ( !i %% 1000 ) {
#    cat( "DBG: Finished ", i, " iterations ", date(), "\n" )
#  }
  entry=loopvar[[i]]
  if( length(entry) > 0 ) {
    fits[[i]] = entry[[1]]

    lambda.i[[i]] = entry[[2]]

    pred[i,] = entry[[3]]

    r.pVa[[i]] = entry[[4]]
    coef.i = entry[[5]]
    match.coef = entry[[6]]
    coefficients[ i, match.coef ] = coef.i
  }
}

#finally, output the TRN model and various perforamce metrics to a list object.
return ( list( fits = fits , predicted.expression = pred , actual.expression = actual , beta.coefficients = coefficients , lambda = lambda.i , r2.predVactual = r.pVa ) )
}
