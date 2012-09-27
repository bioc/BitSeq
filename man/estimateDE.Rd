\name{estimateDE}
\alias{estimateDE}
\title{Estimate condition mean expression and calculate Probability of Positive Log Ratio(PPLR)}
\description{
Estimate condition mean expression for both experimental conditions using the expression estimates obtained by \link{estimateExpression}
}
\usage{
estimateDE(conditions, outFile, parFile, 
      lambda0=NULL, samples=NULL, confidencePerc=NULL, 
      verbose=NULL, norm=NULL, pretend=FALSE )
}
\arguments{
\item{conditions}{\code{List} of vectors, each vector containing names of files containing the expression samples from a replicate (Can be both technical and biological replicates. However, in order to get good results biological replicates for each condition are essential).}
\item{outFile}{Prefix for the output files.}
\item{parFile}{File containing estimated hyperparameters.}
\item{samples}{Produce samples of condition mean expression apart from PPLR and confidence.}
\item{confidencePerc}{Percentage for confidence intervals.}
\item{verbose}{Verbose output.}
Advanced options:
\item{lambda0}{Model parameter lambda_0.}
\item{norm}{Vector of (multiplicative) normalization constants for library size normalization of expression samples. Number of constants has to match the number of expression samples files.} 
\item{pretend}{Do not execute, only print out command line calls for the C++ version of the program.}
}
\details{
This function takes as an input expression samples from biological replicates of two or more conditions and hyperparameters over precision distribution inferred by \code{\link{estimateHyperPar}}.
It uses pseudo-vectors of expression samples from all replicates to infer condition mean expression for each condition. 
The condition mean expression samples are used for computation of the Probability of Positive Log Ratio (PPLR) as well as \eqn{\log_2}{log2} fold change of expression with confidence intervals and average condition mean expression for each transcript. 
Optionally the function can produce also the samples of condition mean expression for each condition.
}
\value{
\item{.pplr}{file containing the PPLR, confidence interval, mean log2 fold change, mean condition mean expressions}
\item{.est}{files containing samples of condition mean expressions for each condition - optional}
\item{.estVar}{file containing samples of inferred variance of the first condition - optional}
}
\author{Peter Glaus}


\seealso{\code{\link{estimateExpression}}, \code{\link{estimateHyperPar}}}
\examples{
\dontrun{
cond1Files = c("data-c0b0.rpkm","data-c0b1.rpkm")
cond2Files = c("data-c1b0.rpkm","data-c1b1.rpkm")
estimateDE(conditions=list(cond1Files, cond2Files), outFile="data.pplr", parFile="data.par", 
           norm=c(1.0, 0.999, 1.0017, 0.9998))
}
}

\keyword{differential expression}
