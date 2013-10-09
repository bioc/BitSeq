\name{estimateVBExpression}
\alias{estimateVBExpression}
\title{Estimate expression of transcripts using VB}
\description{Estimates the expression of transcripts using Variational Bayes inference algorithm}
\usage{
estimateVBExpression (probFile, outFile, outputType=NULL, trInfoFile=NULL, 
      seed=NULL, samples=NULL, optLimit=1e-5, optMethod="FR", procN=4, 
      verbose=FALSE, veryVerbose=FALSE, pretend=FALSE)
}
\arguments{
\item{probFile}{File with alignment probabilities produced by \code{parseAlignment}}
\item{outFile}{Prefix for the output files.}
\item{outputType}{Output type, possible values: \code{theta}, \code{RPKM}, \code{counts}. This is only relevant when the \code{samples} option is used. Default: \code{theta}.}
\item{trInfoFile}{File containing transcript information. (Necessary for RPKM output)}
\item{seed}{Sets the initial random seed for repeatable experiments.}
\item{samples}{Number of samples to be generated from the posterior distribution. Default: no samples are generated.}
\item{verbose}{Verbose output.}
\item{veryVerbose}{Very verbose output.}
\item{procN}{Maximum number of threads to be used. The program will not use more threads that there are MCMC chains.}
Advanced options:
\item{optLimit}{The optimisation limit in terms of minimal gradient or change of bound.}
\item{optMethod}{The optimisation method, use "FR", "HR", or "steepest".}
\item{pretend}{Do not execute, only print out command line calls for the C++ version of the program.}
}
\details{
This function runs Variational Bayes algorithm to estimate the transcript expression. 
The input is the \code{.prob} file containing alignment probabilities which were produced by \code{\link{parseAlignment}}.
Other optional input is the transcript information file specified by \code{trInfoFile} and again produced by \code{parseAlignment}.

It is much faster inference than MCMC which estimates mean expression equally well.
However, the posterior is in form of Dirichlet distribution with underestimated variance.
Use this method in cases when you are only interested in mean expression.
}
\value{
\item{.m_alphas}{file containing mean relative expression of transcripts \eqn{\theta}{theta} and parameters of the Dirichlet distribution. Please note the first line in the file corresponds to the noise transcript.}
If option \code{samples} is used, the program also generates samples based of the \code{outputType}, the default would be file with extension ".VBtheta".
}
\author{Peter Glaus}

\seealso{\code{\link{parseAlignment}}, \code{\link{estimateExpression}}}
\examples{
\dontrun{
setwd(system.file("extdata",package="BitSeq"))
parseAlignment( "data-c0b0.sam", outFile = "data-c0b0.prob", trSeqFile = "ensSelect1.fasta",
      trInfoFile = "data.tr", uniform = TRUE);

estimateVBExpression( probFile="data-c0b0.prob", outFile="data-c0b0-a", outputType="RPKM", 
      samples=1000, trInfoFile="data.tr", seed=47, verbose=TRUE)
estimateVBExpression( probFile="data-c0b0.prob", outFile="data-c0b0-b", trInfoFile="data.tr")
estimateVBExpression( probFile="data-c0b0.prob", outFile="data-c0b0-c", trInfoFile="data.tr", 
      optLimit=1e-6, optMethod = "HS", procN=12, veryVerbose=TRUE);
}
}

\keyword{transcript expression}
