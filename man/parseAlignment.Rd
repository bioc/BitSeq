\name{parseAlignment}
\alias{parseAlignment}
\title{Compute probabilities of alignments}
\description{
  Compute probability of alignments and save them into \emph{.prob} file.
}
\usage{
parseAlignment( alignFile, outFile, trSeqFile, inputFormat=NULL, trInfoFile=NULL, 
      expressionFile=NULL, readsN=NULL, uniform=TRUE, limitA=NULL, lenMu=NULL,
      lenSigma=NULL, verbose=NULL, veryVerbose=NULL, procN=NULL, pretend=FALSE)
}
\arguments{
\item{alignFile}{File containing read alignments.}
\item{outFile}{Name of the output file.}
\item{inputFormat}{Input format: possible values \code{SAM}, \code{BAM}. (This should be detected automatically in most cases.}
\item{trInfoFile}{If transcript reference sequence information is contained within SAM file, program will write this information into <trInfoFile>, otherwise it will look for this information in the <trInfoFile>.}
\item{trSeqFile}{File containing transcript sequence in FASTA format.}
\item{expressionFile}{Transcript relative expression estimates --- for better non-uniform read distribution estimation.}
\item{readsN}{Total number of reads. This is usually not necessary if SAM/BAM contains also reads with no valid alignments.}
%\item{failed}{File name where to save names of reads that failed to align as pair.}
\item{uniform}{Use uniform read distribution.}
\item{limitA}{Limit maximum number of alignments per read. Reads with more alignments than limit will be discarded.}
\item{lenMu}{Set mean of log fragment length distribution. \eqn{l_{frag} \sim LogNormal(\mu,\sigma^2)}{l_frag ~ LogNormal(mu,sigma^2)}}
\item{lenSigma}{Set \eqn{\sigma^2}{sigma^2} (or variance) of log fragment length distribution. \eqn{l_{frag} \sim LogNormal(\mu,\sigma^2)}{l_frag ~ LogNormal(mu,sigma^2)}}
%\item{distributionFile}{Name of file to which read-distribution should be saved.}
\item{verbose}{Verbose output.}
\item{veryVerbose}{Very verbose output.}
\item{procN}{Maximum number of threads to be used.}
\item{pretend}{Do not execute, only print out command line calls for the C++ version of the program.}
}
\details{
This function uses the alignments and reference file to assign probability to each alignment.
It uses either bias-corrected or uniform model for the read distribution, assumes Log-Normal distribution of fragment lengths for pair-end read data and uses quality scores and mismatches to assign probability for every alignment of a read (or fragment) to a transcript.
}
\value{
\item{.prob}{file containing the alignment probabilities}
\item{.tr}{file containing reference transcript names, lengths and effective lengths - optional}
}

\author{Peter Glaus}
\seealso{\code{\link{estimateExpression}}}
\examples{\dontrun{
parseAlignment(alignFile="data.sam", outFile="data.prob", 
      trSeqFile="trReference.fa" ,trInfoFile="data.tr")
}}

\keyword{transcript expression}
