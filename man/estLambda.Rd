\name{estLambda}
\alias{estLambda}
\title{Determine a reasonable glmnet lambda penalty parameter.}
\description{Determine an appropriate lambda parameter by performing a cross-validation with cv.glmnet.}
\usage{estLambda( expr , physnet , n = 100 , lambda = "adaptive"  )}
\arguments{
  \item{expr}{An expression matrix with genes in rows and samples in columns. Gene names should be in the same order as the training data used to generate 'trn'.}
  \item{physnet}{The TFBS matrix used to make trn.}
  \item{n}{The number of genes that should be sampled by estLambda. Default is 100.}
  \item{lambda}{Method for selecting lambda. Ordinarily, this is set to "adaptive" for estLambda.}
}
\value{
  \item{fits}{glmnet fitted objects for each gene. See fitTRN.}
  \item{predicted.expression}{matrix of predicted expression. See fitTRN.}
  \item{actual.expression}{matrix of centered and scaled expression data. See fitTRN.}
  \item{beta.coefficients}{TFs x target genes matrix of beta-coefficients. See fitTRN.}
  \item{lambda}{lambda parameters used for each gene. Determiend by a call to cv.glmnet for each gene. The lambda.1se value from cv.glmnet is used by default.}
  \item{r2.predVactual}{Variance explained by the model (R^2)}
  \item{lambda.median}{The median lambda from all genes tested by estLambda.}
  \item{r2.median}{The median R^2 from all genes tested by estLambda.}
}
\author{Seth Ament}

