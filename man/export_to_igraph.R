\name{export_to_igraph}
\alias{export_to_igraph}
\title{Export TRN edges to an igraph graph object}
\description{The igraph package provides numerous functions for network analysis.}
\usage{export_to_igraph( trn , cols = NULL, rows = NULL , beta.threshold = 0 ) }
\arguments{
  \item{trn }{A transcriptional regulatory network object produced by fitTRN.}
  \item{cols }{a character vector specifying a subset of colnames(trn) }
  \item{rows }{a character vector specifying a subset of rownames(trn) }
  \item{ beta.threshold }{minimum abs( beta coefficient) for selecting edges}
}
\value{An igraph graph object with edges from the selected rows and columns of trn. Produced by a call to the igraph function graph_from_data_frame().}
\author{Seth Ament}
