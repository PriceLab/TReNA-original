show_hubs =
function( trn , n = 10 ) {
        tfs = colnames(trn)
        kOut = colSums( trn != 0 )
        positive_reg = colSums( trn > 0 )
        negative_reg = colSums( trn < 0 )
        combined_stats = data.frame( tfs , kOut , positive_reg , negative_reg )
        tf_order = order( kOut , decreasing = T )[1:n]
        combined_stats[ tf_order , ]
}

