\name{summarystats}
\alias{summarystats}
\title{ Summary Statistics for a TRN Model.}
\description{number of target genes, number of regulators, number of edges, median kIn , median kOut }
\usage{summarystats( trn )}
\arguments{
  \item{trn }{A transcriptional regulatory network object produced by fitTRN(), or a target genes x TFs matrix produced by trimTRN().}
}
\value{
  \item{nTargets }{ the number of target genes in the network }
  \item{nRegulators}{ the number of TFs in the network }
  \item{nEdges }{ the number of links (edges) in the network }
  \item{kIn.median }{ median in-degree per target gene }
  \item{kOut.median }{ median out-degree per regulator }
}
\author{Seth Ament}