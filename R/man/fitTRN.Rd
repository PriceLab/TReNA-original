\name{fitTRN}
\alias{fitTRN}
\title{Fit a Transcriptional Regulatory Network (TRN)}
\description{fitTRN requires two input files: (1) 'physnet', a target genes x transcription factors (TFs) matrix describing binding sites of TFs near potential target genes. (2) a gene expression matrix. Lasso regression (Tibshirani R, 1996. J. R. Statist. Soc. B. 58: 267-288) -- implemented by a call to the glmnet R package -- is used to fit a linear model that predicts the expression of each target gene based on expression levels of TFs.}
\usage{fitTRN(expr, physnet , lambda = 0.1)}
\arguments{
  \item{expr}{A gene expression matrix with a row corresponding to each gene and a column corresponding to each sample. Some or all of the row names in expr must match (some or all) of the column and row names in physnet. Typically, these row names are standard identifiers such as Entrez ID's or official gene symbols. In principle, expression data may be from any source (e.g., log2-normalized expression values from microarrays, FPKM estimates from RNA-seq). fitTRN expects a single row corresponding to each gene. When multiple values per gene are present in the original data (e.g., multiple probes on a microarray), we typically select a single probe or use the average expression of all probes. fitTRN automatically centers and scales the expression data, but other pre-processing should be done prior to running fitTRN.}
  \item{physnet}{A logical (1/0) matrix with a row corresponding to each gene and column corresponding to each transcription factor (TF). Typically, physnet is derived from information about transcription factor binding sites (TFBSs), based on cis-regulatory motifs or other information. Some or all of the column and row names in physnet must match row names of expr. Typically, these row names are standard identifiers such as Entrez ID's or official gene symbols. physnet maps for the human and mouse genome using Entrez IDs and Official Gene Symbols are provided in data().}
  \item{lambda}{The lasso regression penalty, either a numeric or "adaptive". By default, lambda is set to 0.1. See the man page for glmnet. If lambda is numeric, then this value is used for all target genes. If lambda = "adaptive", an appropriate value for lambda is set independently for each target gene, by running an internal cross-validation with cv.glmnet. lambda is then set as the maximum value of lambda such that the deviance is < 1 s.e above the minimum deviance achieved in the full model (see man pacge for cv.glmnet).}
}
\value{
	\item{fits }{A list containing the glmnet fit objects for each row of expr.}
	\item{predicted.expression }{A matrix with the same dimensions as expr containing the predicted values of each gene based on the fitTRN model.}
	\item{actual.expression }{A matrix with the mean-centered and Z-scaled values of expr.}
	\item{beta.coefficients }{A matrix with nrow() = nrow(expr), and with ncol() = ncol(physnet) + 1. Rows are the same as expr. The first column is an intercept (ordinarily close to 0). The remaining columns of coefficients are in the same order as physnet.}
	\item{lambda }{Values of lambda used for each gene.}
	\item{r2.predVactual }{Pearson's r^2, indicating the fraction of variance explained by the fitTRN model for each gene.}
}
\references{Ament S.A. et al. In preparation.}
\author{
Seth Ament, Institute for Systems Biology (2015)
}