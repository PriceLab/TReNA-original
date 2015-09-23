testTRN <-
function( trn , expr , physnet ) {

require( glmnet )

#remove columns of the TFBS (TFs) for which we do not have expression data.
physnet = physnet[ , which( colnames(physnet) %in% rownames(expr) ) ]

fits = trn$fits
n = length(trn$r2.predVactual)
lambda = trn$lambda
actual.train = trn$actual

#center and scale the expression data.
sd = apply( expr , 1 , sd )
actual = ( expr - rowMeans(expr) ) / sd

#initialize a matrix to store the test set predictions
pred = matrix( NA , nrow = nrow( actual ) , ncol = ncol( actual ) )
r.pVa = rep( NA , n )
tf.symbols = intersect( rownames(expr) , colnames( physnet ) )
pred = matrix( NA , nrow = nrow( actual ) , ncol = ncol( actual ) )

for( i in 1:n ) {

if( i > length(fits) ) next

fit = fits[[i]]
s = lambda[i]

#symbol is the ID of the target gene.
symbol = rownames(actual)[i]

#print the name of the target gene on the console
status = paste( "Working on row" , i , "," , symbol )
cat(status)

#y is the expression data for the target gene
y = actual[ i , ]

#sub defines the row for this target gene in the TFBS database
sub = which( rownames(physnet) == symbol )

#if the target gene is not present in the TFBS database, we skip this gene.
if( length(sub) == 0 ) {
cat("\n")
next
}

#most genes are represented by a single row
#has_tfbs defines the columns (TFs) that have a binding site in the promoter of this target gene.
if( length(sub) == 1 ) {
has_tfbs = which( physnet[sub,] > 0 )
potential.regulators = names(has_tfbs)
}

#a few edge cases may have multiple rows in the TFBS database. We take the union of TFs with binding sites.
if( length(sub) > 1 ) {
has_tfbs = which( colSums(physnet[sub,]) > 0 )
potential.regulators = names(has_tfbs)
}

#potential.regulators is a vector with the names of TFs that have binding sites in the promoter of target gene i#
#we consider only those TFs for which we have expression data, defined by the vector tf.symbols.
potential.regulators = intersect( potential.regulators , tf.symbols )

#if there are no predicted regulators, we skip this gene.
if( length(potential.regulators) == 0 ) cat("\n")
if( length(potential.regulators) == 0 ) next

#autoregulation -- in which a TF regulates itself -- throws off the model. Although autoregulation is potentially interesting, we eliminate autoregulatory loops.
pr = setdiff( which( rownames(actual) %in% potential.regulators) , which( rownames(actual)==symbol) )

#if there are no predicted regulators after eliminating autoregulatory loops, we discard this gene.
if( length(pr) < 2 ) cat("\n")
if( length(pr) < 2 ) next

#x is a matrix with the expression levels of TFs (potential regulators) in each of the samples.
x = t( actual[ pr , ] )

pred.i = predict( fit , newx = x , s = s , type = "link" )
pred[i,] = pred.i
r.pVa[i] = cor( pred.i[,1] , y )

cat( "\n" )

}

list( predicted.expression = pred , r2.predVactual = r.pVa ^ 2 )

}
