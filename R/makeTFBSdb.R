makeTFBSdb = function( 	x , 
			motif_to_tf , 
			p.fimo.thresh = 1 , 
			p.wellington.thresh = 1 ,
			overlap.thresh = 1 , 
			dist_to_tss.max = 10000 ,
			dist_to_tss.min = 0 ,
			verbose = T ,
			id = "symbol" ) {

if( verbose == T ) cat( "filtering motifs\n" )
high_quality = which(
	x$p.fimo < p.fimo.thresh &
	x$p.wellington < p.wellington.thresh &
	x$motif_footprint_overlap_bp >= overlap.thresh &
	abs(x$dist_to_tss) < dist_to_tss.max & 
	abs(x$dist_to_tss) >= dist_to_tss.min )
motif_table = x[ high_quality , ]

if( verbose == T ) cat( "calculating target genes for each TF\n" )
tfs = unique( motif_to_tf[,2] )

if( id == "symbol" ) {
 genes = unique( motif_table$symbol )
}
if( id == "ensembl_gene_id" ) {
 genes = unique( motif_table$ensembl_gene_id )
}

ntf = length(tfs)
ngene = length(genes)
# main loop
tfbsdb = matrix( 0 , nr=ngene , nc=ntf )
colnames(tfbsdb) = tfs
rownames(tfbsdb) = genes
for( i in 1:ntf ) {
 motifs.i = unique( motif_to_tf[ which( motif_to_tf[,2] == tfs[i] ) , 1 ])
 instances = which( motif_table$motif %in% motifs.i )
 gene_counts = table( motif_table[ instances , "symbol" ] )
 matchgenes = match( names(gene_counts) , genes )
 tfbsdb[ matchgenes , i ] = gene_counts
 # cat( i , "\n" )
} # end main loop

# summary statistics
kin = rowSums( tfbsdb > 0 )
kout = colSums( tfbsdb > 0 )
q.kin = quantile( kin , c(.25,.5,.75) )
q.kout = quantile( kout , c(.25,.5,.75) )
if( verbose == T ) {
  cat(  "done" , i , "tfs\n" ,
	"summary statistics:\n" ,
	"\t0%\t25%\t50%\t75%\t100%\n" ,
	"kIn\t" , min(kin) , "\t" , q.kin[1],"\t" , q.kin[2],"\t" , q.kin[3],"\t" , max(kin) , "\n" ,
        "kOut\t" , min(kout) , "\t" , q.kout[1],"\t" , q.kout[2],"\t" , q.kout[3],"\t" , max(kout) , "\n" )
}

return( tfbsdb )

}

