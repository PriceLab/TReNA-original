\name{setEnrich}
\alias{setEnrich}
\title{Calculate enrichments of TRN modules for gene sets.}
\description{Fisher's exact test enrichments of the targets of each TF for the genes in a user-defined set (e.g., differentially-expressed genes, a GO term)}
\usage{setEnrich trn , set , minModSize = 15 , universe = NULL , verbose = T )}
\arguments{
  \item{trn }{A transcriptional regulatory network object produced by fitTRN.}
  \item{set }{The gene set. A character vector corresponding row names in trn (e.g., gene symbol).}
  \item{minModSize }{minimum size of TRN modules to be tested}
  \item{ universe }{character vector giving the background ("universe") of genes for 'set'. If universe = NULL (the default), rownames(trn) is used as the background for enrichment calculations.}
  \item{ verbose }{ Logical. Should information about the progress of calculations be printed to the console? }
}
\value{A data frame with a row for each column (TF) in trn, and the following columns: tf, the TF name; overlap, the number of genes shared between the TF's targets and set; expect, the size of overlap expected by chance; rf, relative frequency (overlap / expect); p, Fisher's exact test p-value; q, Benjamini-Hochberg FDR across all TFs.}
\author{Seth Ament}
