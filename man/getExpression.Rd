\name{getExpression}
\alias{getExpression}
\title{Estimate transcript expression}
\description{
Estimate expression of transcripts. Starting from alignment and reference files function function handles the entire process of expression analysis resulting in transcript expression means and variances together with file containing all the expression samples.
}
\usage{
getExpression(alignFile, trSeqFile, outPrefix=NULL, uniform=TRUE, type="RPKM", log=FALSE, ... )
}
\arguments{
\item{alignFile}{File containing read alignments.}
\item{trSeqFile}{File containing transcript sequence in FASTA format.}
\item{outPrefix}{Prefix for the output files. Otherwise program creates temporary files, which are only valid for current R session.}
\item{uniform}{Use uniform read distribution.}
\item{type}{Output type, possible values: \code{theta}, \code{RPKM}, \code{counts}, \code{tau}.}
\item{log}{Report mean and expression of logged expression samples.}
\item{...}{Other arguments are passed to \code{\link{estimateExpression}}, for more details see \code{\link{estimateExpression}}}
}
\details{
This function uses \code{parseAlignment} function to compute alignment probabilities and the function \code{estimateExpression} to produce the expression samples.

In case of non-uniform read distribution, it first produces approximate estimates of expression using uniform distribution and uses these estimates in to compute read distribution bias-corrected alignment probabilities, which are used in the \code{estimateExpression} function to produce expression estimates.
}
\value{
\code{list} with items:
\item{exp}{\code{DataFrame} with transcript expression mean and variance}
\item{fn}{name of the file containing all the expression samples}
}
\author{Peter Glaus}

\seealso{\code{\link{getDE}}, \code{\link{estimateExpression}}, \code{\link{parseAlignment}}}
\examples{\dontrun{
res1 <- getExpression("data-c0b0.sam","ensSelect1.fasta", MCMC_chains=2, MCMC_samplesN=100)
}}

\keyword{transcript expression}
