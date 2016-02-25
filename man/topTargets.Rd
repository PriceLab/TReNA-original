\name{topTargets}
\alias{topTargets}
\title{for a given TF, print the top n predicted targets}
\description{for a given TF, print the top n predicted targets}
\usage{topTFs(trn,tf,n=10)}
\arguments{
  \item{trn}{A TRN object, usually produced by fitTRN.}
  \item{gene}{either a character string corresponding to the name of a TF (a column of trn), or an integer corresponding to a column of trn.}
  \item{n}{the number of target genes to show. By default, show the top 10 target genes, ordered by the strength of the beta coefficients.}
}
\value{Target gene names and coefficients are printed to the console.}
\author{Seth Ament}
