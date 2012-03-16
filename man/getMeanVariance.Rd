\name{getMeanVariance}
\alias{getMeanVariance}
\title{Calculate mean and variance of expression samples}
\description{Calculate mean and variance of expression samples or log-expression samples}
\usage{
getMeanVariance(sampleFiles, outFile, log=NULL, type=NULL, verbose=NULL)
}
\arguments{
\item{sampleFiles}{List of one or more files containing the expression samples.}
\item{outFile}{Name of the output file.}
\item{log}{Use logged values.}
\item{type}{Type of variance, possible values: \code{sample},\code{sqDif} for sample variance or squared difference.}
\item{verbose}{Verbose output.}
}
\details{
The \code{getMeanVariance} function computes means and variances of MCMC expression samples.
These can be computed either from single file or from multiple files using sample variance. 
Variance of two experiments (i.e. technical or biological replicates) can be estimated also by using \code{sqDif} option for \code{type} which specify the computation of the average square distance between the samples from two sets.
}
\value{
\item{.means}{File containing means (first column) and variance (second column) for each transcript (or row in the sample files)}
}
\author{Peter Glaus}

\seealso{\code{\link{estimateExpression}}}
\examples{\dontrun{
sampleFileNames = c("data-c0b0.rpkm","data-c0b1.rpkm","data-c1b0.rpkm","data-c1b1.rpkm")
getMeanVariance(sampleFiles=sampleFileNames, outFile="data.means", log=1)
}}

\keyword{expression mean}
