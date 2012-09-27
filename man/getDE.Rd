\name{getDE}
\alias{getDE}
\title{Estimate Probability of Positive Log Ratio}
\description{Using expression samples, program estimates the probability of differential expression for each transcript.}
\usage{
getDE(conditions, outPrefix=NULL, samples=FALSE, trInfoFile=NULL, 
      norm=NULL, pretend=FALSE )
}
\arguments{
\item{conditions}{\code{List} of vectors, each vector containing names of files containing the expression samples from a replicate (Can be both technical and biological replicates. However, in order to get good results biological replicates for each condition are essential).}
\item{outPrefix}{Prefix for the output files. Otherwise program creates temporary files, which are only valid for current R session.}
\item{samples}{Produce samples of condition mean expression apart from PPLR and confidence.}
\item{trInfoFile}{Transcript information file providing the names of transcripts.}
\item{norm}{Vector of (multiplicative) normalization constants for library size normalization of expression samples. Number of constants has to match the number of expression samples files.} 
\item{pretend}{Do not execute, only print out command line calls for the C++ version of the program.}
}
\details{
This function uses \code{estimateHyperPar} function to estimate the hyperparameters for DE model and the uses \code{estimateDE} function to infer the condition mean expression and calculate Probability of Positive Log Ratio.
}
\value{
   \code{list} with items:
   \item{pplr}{\code{DataFrame} with PPLR and other statistics}
   \item{fn}{\code{list} with file names for PPLR file, \code{fn$pplr}, and condition mean expression samples, \code{fn$samplesFiles} (only with option \code{samples=TRUE})}
}
\author{Peter Glaus}

\seealso{\code{\link{getExpression}}, \code{\link{estimateHyperPar}}, \code{\link{estimateDE}}}
\examples{\dontrun{
cond1Files = c("data-c0b0.rpkm","data-c0b1.rpkm")
cond2Files = c("data-c1b0.rpkm","data-c1b1.rpkm")
deRes <- getDE( conditions=list(cond1Files, cond2Files))
## top 10 DE transcripts
head(deRes$pplr[ order(abs(0.5-deRes$pplr$pplr), decreasing=TRUE ), ], 10)
}}

\keyword{differential expression}
