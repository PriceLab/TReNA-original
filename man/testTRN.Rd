\name{testTRN}
\alias{testTRN}
\title{Test accuracy of predictions from fitTRN in independent samples.}
\description{Uses the lasso regression model from a TRN object (usually produced by fitTRN) to predict expression levels in an independent test set.}
\usage{testTRN(trn, expr, physnet)}
\arguments{
  \item{trn}{A TRN object, usually produced by fitTRN.}
  \item{expr}{An expression matrix with genes in rows and samples in columns. Gene names should be in the same order as the training data used to generate 'trn'.}
  \item{physnet}{The TFBS matrix used to make trn.}
}
\value{
	\item{predicted.expression }{A matrix with the same dimensions as expr containing the predicted values of each gene based on the fitTRN model.}
	\item{r2.predVactual }{Pearson's r^2, indicating the fraction of variance explained by the fitTRN model for each gene.}
}
\author{Seth Ament}
