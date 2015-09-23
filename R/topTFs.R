topTFs =
function( trn , gene , n = Inf ) {
        if( is.integer( gene ) ) gene_index = gene
        if( is.character( gene ) ) gene_index = which( rownames( trn ) == gene )
        gene_name = rownames( trn ) [ gene_index ]
        regulator_indices = which( trn[ gene_index , ] != 0 )
        regulator_coefs = trn[ gene_index , regulator_indices ]
        abs_coefs = abs( regulator_coefs )
        regulator_names = colnames( trn[ , regulator_indices ] )
        kIn = length( regulator_indices )
        regulator_rank = order( abs_coefs , decreasing = T )
        n_used = min( kIn , n )
        cat( paste( "Top Regulators for" , gene_name , "\n" ))
        cat( paste( "Showing top" , n_used , "of" , kIn , "regulators.\n" ))
        cat( paste( "   " , "Symbol" , "\t" , "TRN Coefficient" , "\n" ))
        for( i in regulator_rank[1:n_used] ) {
                line = paste( "   " , regulator_names[i] , "\t" , round(regulator_coefs[i],2) , "\n" )
                cat(  line )
        }
}

