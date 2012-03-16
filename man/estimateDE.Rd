\name{estimateDE}
\alias{estimateDE}
\title{Estimate condition mean expression and calculate Probability of Positive Log Ratio(PPLR)}
\description{
Estimate condition mean expression for both experimental conditions using the expression estimates obtained by \link{estimateExpression}
}
\usage{
estimateDE( cond1, cond2, outFile, parFile, lambda0=NULL, samples=NULL, confidencePerc=NULL, verbose=NULL )
}
\arguments{
\item{cond1}{List of files containing the expression samples for contition 1.}
\item{cond2}{List of files containing the expression samples for contition 2.}
\item{outFile}{Prefix for the output files.}
\item{parFile}{File containing estimated hyperparameters.}
\item{samples}{Produce samples of condition mean expression apart from PPLR and confidence.}
\item{confidencePerc}{Percentage for confidence intervals.}
\item{verbose}{Verbose output.}
Advanced options:
\item{lambda0}{Model parameter lambda_0.}
}
\details{
This function takes as an input expression samples from biological replicates of two conditions and hyperparameters over precision distribution inferred by \code{\link{estimateHyperPar}}.
It uses pseudo-vectors of expression samples from all replicates to infer condition mean expression for each condition. 
The condition mean expression samples are used for computation of the Probability of Positive Log Ratio (PPLR), with confidence intervals as well as \eqn{\log_2}{log2} fold change of expression and average condition mean expression for each transcript. 
Optionally the function can produce also the samples of condition mean expression for each condition.
}
\value{
\item{.pplr}{file containing the PPLR, confidence interval, mean log2 fold change, mean condition mean expressions}
\item{.est}{files containing samples of condition mean expressions for each condition - optional}
\item{.estVar}{file containing samples of inferred variance of the first condition - optional}
}
\author{Peter Glaus}


\seealso{\code{\link{estimateExpression}}, \code{\link{estimateHyperPar}}}
\examples{\dontrun{
cond1Files = c("data-c0b0.rpkm","data-c0b1.rpkm")
cond2Files = c("data-c1b0.rpkm","data-c1b1.rpkm")
estimateDE( cond1=cond1Files, cond2=cond2Files, outFile="data.pplr", parFile="data.par" )
}}

\keyword{differential expression}
