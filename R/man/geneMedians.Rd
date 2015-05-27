\name{geneMedians}
\alias{geneMedians}
\title{Generate a genes x samples expression matrix from probe- or transcript-level expression data.}
\description{Generate a genes x samples expression matrix from probe- or transcript-level expression data by taking the median expression of all the probes for each gene.}
\usage{geneMedians( expr, anno )}
\arguments{
  \item{expr}{An expression matrix with features in rows and samples in columns.}
  \item{anno}{A character vector of length nrow(expr) mapping rows of expr to the names of genes.}
}
\value{A matrix with a row for each unique gene identifier in anno and columns matching the samples in expr.}
\author{Seth Ament}