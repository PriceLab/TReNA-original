\name{trimTRN}
\alias{trimTRN}
\title{Trim a transcriptional regulatory network model.}
\description{Remove poorly predicted genes, inactive TFs, and weak connections from a TRN model produced by fitTRN, using user-defined thresholds.}
\usage{trimTRN( trn , r2.threshold = 0.5 , beta.threshold = 0.01 )}
\arguments{
  \item{trn }{A transcriptional regulatory network object produced by fitTRN.}
  \item{r2.threshold }{Inclusion threshold for prediction accuracy of each target gene, i.e. the fraction of expression variance explained by the model. Default is 0.5.}
  \item{beta.threshold }{Threshold below which edges with small predicted effect sizes should be excluded from the model. Default is 0.01.}
}
\value{A numeric matrix with a row for each accurately-predicted gene, a column for each TF with >0 predicted target genes, and beta coefficients for the predicted effect of the TFs on each target gene.}
\author{Seth Ament}