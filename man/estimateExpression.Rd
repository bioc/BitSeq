\name{estimateExpression}
\alias{estimateExpression}
\alias{estimateExpressionLegacy}
\title{Estimate expression of transcripts}
\description{Estimates the expression of transcripts using Markov chain Monte Carlo Algorithm}
\usage{
estimateExpression(probFile, outFile, parFile=NULL, outputType=NULL, gibbs=NULL, 
      trInfoFile=NULL, thetaActFile=NULL, MCMC_burnIn=NULL, MCMC_samplesN=NULL, 
      MCMC_samplesSave=NULL, MCMC_chainsN=NULL, MCMC_dirAlpha=NULL, seed=NULL, 
      verbose=NULL, procN=NULL, pretend=FALSE)
estimateExpressionLegacy(probFile, outFile, parFile=NULL, outputType=NULL, gibbs=NULL, 
      trInfoFile=NULL, thetaActFile=NULL, MCMC_burnIn=NULL, MCMC_samplesN=NULL, 
      MCMC_samplesSave=NULL, MCMC_samplesNmax=NULL, MCMC_chainsN=NULL, 
      MCMC_scaleReduction=NULL, MCMC_dirAlpha=NULL, seed=NULL, verbose=NULL, 
      pretend=FALSE)
}
\arguments{
\item{probFile}{File with alignment probabilities produced by \code{parseAlignment}}
\item{outFile}{Prefix for the output files.}
\item{outputType}{Output type, possible values: \code{theta}, \code{RPKM}, \code{counts}, \code{tau}.}
\item{gibbs}{Use regular Gibbs sampling instead of Collapsed Gibbs sampling.}
\item{parFile}{File containing parameters for the sampler, which can be otherwise specified by [MCMC*] options. As the file is checked after every MCMC iteration, the parameters can be adjusted while running.}
\item{trInfoFile}{File containing transcript information. (Necessary for RPKM)}
\item{MCMC_burnIn}{Length of sampler's burn in period.}
\item{MCMC_samplesN}{Initial number of samples produced. These are used either to estimate the number of necessary samples or to estimate possible scale reduction.}
\item{MCMC_samplesSave}{Number of samples recorder at the end in total.}
\item{MCMC_chainsN}{Number of parallel chains used. At least two chains will be used.}
\item{seed}{Sets the initial random seed for repeatable experiments.}
\item{verbose}{Verbose output.}
\item{procN}{Maximum number of threads to be used. The program will not use more threads that there are MCMC chains.}
Advanced options:
\item{thetaActFile}{File for logging noise parameter thetaAct, which is only generated when regular Gibbs sampling is used.}
\item{MCMC_dirAlpha}{Alpha parameter for the Dirichlet distribution.}
\item{pretend}{Do not execute, only print out command line calls for the C++ version of the program.}
\item{MCMC_scaleReduction}{(Only for \code{estimateExpressionLegacy}.) Target scale reduction, sampler finishes after this value is met.}
\item{MCMC_samplesNmax}{(Only for \code{estimateExpressionLegacy}.) Maximum number of samples produced in one iteration. After producing samplesNmax samples sampler finishes.}
}
\details{
This function runs Collapse Gibbs algorithm to sample the MCMC samples of transcript expression. 
The input is the \code{.prob} file containing alignment probabilities which were produced by \code{\link{parseAlignment}}.
Other optional input is the transcript information file specified by \code{trInfoFile} and again produced by \code{parseAlignment}.

The \code{estimateExpression} function first runs burn-in phase and initial iterations to estimate the properties of the MCMC sampling. 
The initial samples are used to estimate the number of samples necessary for generating \code{MCMC_samplesSave} effective samples in the second, final, stage.

The \code{estimateExpressionLegacy} uses less efficient convergence checking via "scale reduction" estimation. 
After an iteration of generating \code{MCMC_samplesN} samples, it estimates possible scale reduction of the marginal posterior variance. 
While the possible scale reduction is high, it doubles the \code{MCMC_samplesN} and starts new iteration. 
This process is repeated until desired value of \code{MCMC_scaleReduction} is met, or \code{MCMC_samplesNmax} samples are generated.

The sampling algorithm can be configured via parameters file \code{parFile} or by using the \code{MCMC*} options.
The advantage of using the file (at least an existing blank text document) is that by changing the configuration values while running, the new values do get updated after every iteration. 
}
\value{
\item{.thetaMeans}{file containing average relative expression of transcripts \eqn{\theta}{theta}}
Either one of sample files based on output type selected:
\item{.rpkm}{ for RPKM expression}
\item{.counts}{ for estimated read counts}
\item{.theta}{ for relative expression of fragments}
\item{.tau}{ for relative expression of transcripts}
}
\author{Peter Glaus}

\seealso{\code{\link{parseAlignment}}}
\examples{\dontrun{
estimateExpression( probFile="data.prob", outFile="data", outputType="RPKM", 
      trInfoFile="data.tr", seed=47, verbose=TRUE)
estimateExpression( probFile="data-c0b0.prob", outFile="data-c0b0", outputType="RPKM", 
      trInfoFile="data.tr", MCMC_burnIn=200, MCMC_samplesN=200, MCMC_samplesSave=100, 
      MCMC_chainsN=2 , MCMC_dirAlpha=NULL )
estimateExpression( probFile="data.prob", outFile="data-G", gibbs=TRUE, 
      parFile="parameters1.txt", outputType="counts", trInfoFile="data.tr")
estimateExpressionLegacy( probFile="data-c0b0.prob", outFile="data-c0b0", outputType="RPKM", 
      trInfoFile="data.tr", MCMC_burnIn=200, MCMC_samplesN=200, MCMC_samplesSave=100, 
      MCMC_samplesNmax=10000, MCMC_scaleReduction=1.2, MCMC_chainsN=2 , MCMC_dirAlpha=NULL )
}}

\keyword{transcript expression}
