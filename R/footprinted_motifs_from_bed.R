footprinted_motifs_from_bed = function( bed , windowsize = 1000000 ) {

# load the bedfile with overlap between output from 
# FIMO, Wellington, and ENSEMBL ranges +-100kb from TSSs
colnames(bed) = c(
	"chrom.gene" ,
	"tss_minus_window" ,
	"tss_plus_window" ,
	"ensembl_gene_id" ,
	"symbol" ,
	"chrom.motif" ,
	"start.motif" ,
	"end.motif" ,
	"score.motif" , 
	"strand.motif" ,
	"motif" ,
	"motif.description1" ,
	"motif.description2" ,
	"unnamed_overlap1" ,
	"chrom.footprint" ,
	"start.footprint" ,
	"end.footprint" ,
	"footprint_name" ,
	"logp_footprint" ,
	"strand.footprint" ,
	"motif_footprint_overlap_bp" ,
	"unnamed_overlap2" )

# calculate the distance from the motif instance to the transcription start site
tss = bed$tss_plus_window - windowsize
dist_to_tss = bed$start.motif - tss

# extract the motif names
# motif = gsub( "^(.*?)\\=" , "" , bed$motif.description1 )
# motif = gsub( "\\;(.*)" , "" , motif )

motif = bed$motif

# extract the motif p-value
p.fimo = gsub( "^(.*?)pvalue=" , "" , bed$motif.description1 )
p.fimo = gsub( "\\;(.*)" , "" , p.fimo )
p.fimo = as.numeric( p.fimo )

# extract footprint p-value
p.wellington = 10^(bed$logp_footprint)
 
# create a new object with the most important information
footprinted_motifs = 
	data.frame( 	
	chrom = bed$chrom.motif ,
	start = bed$start.motif ,
	end = bed$end.motif ,
	motif = motif ,
	ensembl_gene_id = bed$ensembl_gene_id ,
	symbol = bed$symbol ,
	dist_to_tss = dist_to_tss ,
	p.fimo = p.fimo ,
	p.wellington = p.wellington ,
	motif_footprint_overlap_bp = bed$motif_footprint_overlap_bp 
	)

return( footprinted_motifs )

}

