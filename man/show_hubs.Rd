\name{show_hubs}
\alias{show_hubs}
\title{Print the most connected TFs in a TRN model.}
\description{Compute kOut, number of positively regulated, and number of negatively regulated tartgets for each TF.}
\usage{show_hubs( trn , n = 10 )}
\arguments{
  \item{trn}{A TRN object, usually produced by fitTRN.}
  \item{n}{Number of TFs to show in output}
}
\value{ A data frame with four columns: TF name, kOut , number of positively regulated targets , number of negatively regulated targets.}
\author{Seth Ament}
