\name{estimateHyperPar}
\alias{estimateHyperPar}
\title{Estimate hyperparameters for DE model using expression samples and joint mean expression}
\description{Estimate hyperparameters for the Differential Expression model using expression samples and produced smoothed values of the hyperparameters depending on joint mean expression.}
\usage{
estimateHyperPar( outFile, conditions=NULL, paramsInFile=NULL, 
      meanFile=NULL, force=TRUE, exThreshold=NULL, lambda0=NULL, 
      paramsAllFile=NULL, smoothOnly=NULL, lowess_f=NULL, lowess_steps=NULL, verbose=NULL, 
      veryVerbose=NULL, norm=NULL, pretend=FALSE )
}
\arguments{
\item{outFile}{Name of the output file.}
\item{conditions}{\code{List} of vectors, each vector containing names of files containing the expression samples from a replicate (Can be both technical and biological replicates. However, in order to get good results biological replicates for each condition are essential).}
\item{paramsInFile}{File produced by previous run of the function using \code{paramsAllFile} flag.}
\item{meanFile}{Name of the file containing joint mean and variance.}
\item{exThreshold}{Threshold of lowest expression for which the estimation is done.}
\item{paramsAllFile}{Name of the file to which to store all parameter values generated prior to lowess smoothing(good for later, more careful re-smoothing.)}
\item{smoothOnly}{Input file contains previously sampled hyperparameters which should smoothed only.}
\item{verbose}{Verbose output.}
Advanced options:
\item{force}{Force smoothing hyperparameters, otherwise program might not produce parameters file at the end.}
\item{lambda0}{Model parameter lambda0.}
\item{lowess\_f}{Parameter F for lowess smoothing specifying amount of smoothing.}
\item{lowess\_steps}{Parameter Nsteps for lowess smoothing specifying number of iterations.}
\item{veryVerbose}{More verbose output.}
\item{norm}{Vector of (multiplicative) normalization constants for library size normalization of expression samples. Number of constants has to match the number of expression samples files.} 
\item{pretend}{Do not execute, only print out command line calls for the C++ version of the program.}
%\item{groupsNumber}{Number of groups of transcript of similar size.}
%\item{samplesNumber}{Number of samples generated for each group.}
}
%\details{}
\value{
\item{.par}{file containing the smoothed hyperparameters}
\item{.ALLpar}{file containing all hyperparameter samples prior to smoothing - optional}
}
\author{Peter Glaus}

\seealso{\code{\link{estimateDE}}}
\examples{\dontrun{
cond1Files = c("data-c0b0.rpkm","data-c0b1.rpkm")
cond2Files = c("data-c1b0.rpkm","data-c1b1.rpkm")
estimateHyperPar( conditions=list(cond1Files, cond2Files), outFile="data.par", 
      meanFile="data.means", verbose=TRUE)

estimateHyperPar( conditions=list(cond1Files, cond2Files), outFile="data.par", 
      meanFile="data.means", paramsFile="data.ALLpar", force=FALSE)
estimateHyperPar( outFile="data.par", paramsInFile="data.ALLpar", smoothOnly=TRUE )
}}

\keyword{DE model hyperparameters}
