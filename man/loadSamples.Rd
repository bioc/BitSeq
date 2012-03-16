\name{loadSamples}
\alias{loadSamples}
\alias{writeSamples}
\title{Loading and saving expression samples}
\description{Functions for loading expression samples into DataFrame and saving samples from DataFrame into a file.}
\usage{
loadSamples(fileName, trInfoFile=NULL)
writeSamples(data, fileName)
}
\arguments{
\item{fileName}{Name of the file with samples or to which the samples are written.}
\item{data}{DataFrame with samples written to the file.}
\item{trInfoFile}{Transcript information file which can be used to name the rows.}
}
\details{
The \code{loadSamples} function load samples from the specified file into a \code{DataFrame}.
If the transcript information file is provided, the transcript names are use as row names.

The \code{writeSamples} function can save samples from a \code{DataFrame} into a file in format which is valid for BitSeq and can be used in other functions.
}
\value{
\item{DataFrame}{Containing the expression samples}
}
\author{Peter Glaus}

\seealso{\code{\link{estimateExpression}}}
\examples{\dontrun{
samples1<-loadSamples("data-c0b1.rpkm")
writeSamples(samples1,"new-c0b1.rpkm")
}}

\keyword{expression samples}
