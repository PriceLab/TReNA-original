topTargets =
function( trn , tf , n = 10 ) {
        if( is.integer( tf ) ) tf_index = tf
        if( is.character( tf ) ) tf_index = which( colnames( trn ) == tf )
        regulator_name = colnames( trn ) [ tf_index ]
        target_indices = which( trn[ , tf_index ] != 0 )
        target_coefs = trn[ target_indices , tf_index ]
        abs_coefs = abs( target_coefs )
        target_names = rownames( trn[ target_indices , ] )
        kOut = length( target_indices )
        target_rank = order( abs_coefs , decreasing = T )
        n_used = min( kOut , n )
        cat( paste( "Top Targets for" , regulator_name , "\n" ))
        cat( paste( "Showing top" , n_used , "of" , kOut , "target genes.\n" ))
        cat( paste( "   " , "Symbol" , "\t" , "TRN Coefficient" , "\n" ))
        for( i in target_rank[1:n_used] ) {
                line = paste( "   " , target_names[i] , "\t" , round(target_coefs[i],2) , "\n" )
                cat(  line )
        }
}

