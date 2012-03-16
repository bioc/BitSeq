\name{getDE}
\alias{getDE}
\title{Estimate Probability of Positive Log Ratio}
\description{Using expression samples, program estimates the probability of differential expression for each transcript.}
\usage{
getDE( cond1, cond2, outPrefix=NULL, samples=FALSE, trInfoFile=NULL )
}
\arguments{
\item{cond1}{List of files containing the expression samples for condition 1.}
\item{cond2}{List of files containing the expression samples for condition 2.}
\item{outPrefix}{Prefix for the output files. Otherwise program creates temporary files, which are only valid for current R session.}
\item{samples}{Produce samples of condition mean expression apart from PPLR and confidence.}
\item{trInfoFile}{Transcript information file providing the names of transcripts.}
}
\details{
This function uses \code{estimateHyperPar} function to estimate the hyperparameters for DE model and the uses \code{estimateDE} function to infer the condition mean expression and calculate Probability of Positive Log Ratio.
}
\value{
   \code{list} with items:
   \item{pplr}{\code{DataFrame} with PPLR and other statistics}
   \item{fn}{\code{list} with file names for PPLR file \code{pplr} and condition mean expression samples \code{C1samples}, \code{C2samples} (only with option \code{samples=TRUE})}
}
\author{Peter Glaus}

\seealso{\code{\link{getExpression}}, \code{\link{estimateHyperPar}}, \code{\link{estimateDE}}}
\examples{\dontrun{
cond1Files = c("data-c0b0.rpkm","data-c0b1.rpkm")
cond2Files = c("data-c1b0.rpkm","data-c1b1.rpkm")
deRes <- getDE( cond1=cond1Files, cond2=cond2Files )
## top 10 DE transcripts
head(deRes$pplr[ order(abs(0.5-deRes$pplr$pplr), decreasing=TRUE ), ], 10)
}}

\keyword{differential expression}
