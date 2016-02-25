\name{topTFs}
\alias{topTFs}
\title{for a given target gene, print the top n regulators}
\description{for a given target gene, print the top n regulators}
\usage{topTFs(trn,gene,n=Inf)}
\arguments{
  \item{trn}{A TRN object, usually produced by fitTRN.}
  \item{gene}{either a character string corresponding to the name of a gene (a row of trn), or an integer corresponding to a row of trn.}
  \item{n}{the number of TFs to show. By default, show all TFs with non-zero coefficients.}
}
\value{Regulator names and coefficients are printed to the console.}
\author{Seth Ament}
