\name{getGeneExpression}
\alias{getGeneExpression}
\alias{getWithinGeneExpression}
\title{Calculate gene expression or relative within gene expression}
\description{Calculate either gene expression or relative within gene expression using transcript expression samples and transcript information file.}
\usage{
getGeneExpression(sampleFile, outFile=NULL, trInfo=NULL, trInfoFile=NULL,
      pretend=FALSE)
getWithinGeneExpression(sampleFile, outFile=NULL, trInfo=NULL, trInfoFile=NULL,
      pretend=FALSE, keepOrder=FALSE)
}
\arguments{
\item{sampleFile}{File containing the transcript expression samples.}
\item{outFile}{Name of the output file. If not used, function uses temporary file.}
\item{trInfo}{\code{DataFrame} containing transcript information. Either \code{trInfo} or \code{trInfoFile} argument has to be provided. Otherwise function tries file with same name as \code{sampleFile} and extension \code{tr}. }
\item{trInfoFile}{Transcript information file. Either \code{trInfo} or \code{trInfoFile} argument has to be provided. Otherwise function tries file with same name as \code{sampleFile} and extension \code{tr}.}
\item{pretend}{Do not execute, only print out command line calls for the C++ version of the program.}
\item{keepOrder}{If \code{TRUE} then transcripts will always keep same order, otherwise transcripts might be grouped by genes in the output. (The order is always same if transcripts are grouped by genes.)}
}
\details{
The \code{getGeneExpression} function takes samples of transcript expression and produces file with expression of genes by adding up transcript expression.

The \code{getWithinGeneExpression} function takes samples of transcript expression and produces file with relative within gene expression samples for each transcript.

Both function need valid transcript information which contains gene transcript mapping.
This can be provided either via \code{DataFrame} \code{trInfo} or file named \code{trInfoFile}.

In case of a file, it should be formatted in following manner.
The first line should contain "# M <numberOfTranscripts>" and the following \code{numberOfTranscripts} lines have to contain "<geneName> <transcriptName> <transcriptLength>". 
Example is provided in \code{extdata/ensSelect1.tr}. 
Please note that the transcript information file automatically generated from alignment files are not sufficient because SAM/BAM files do not include gene names. 
We hope to provide more convenient way in future versions of BitSeq.
}
\value{
Name of file containing the new expression samples.
}
\author{Peter Glaus}

\seealso{\code{\link{getExpression}}, \code{\link{tri.load}}, \code{\link{tri.file.setGeneNames}}, \code{\link{tri.file.hasGeneNames}}}
\examples{
setwd(system.file("extdata",package="BitSeq"))
## use transcript information as object
trinfo <- tri.load("ensSelect1.tr")
## gene expression
getGeneExpression("data-c0b1.rpkm", "data-c0b1-GE.rpkm", trInfo=trinfo)
gExpSamples <- loadSamples("data-c0b1-GE.rpkm")
gExpMeans <- rowMeans(as.data.frame(gExpSamples))
gExpMeans

## within gene expression
wgeFN <- getWithinGeneExpression("data-c0b1.rpkm", trInfoFile="ensSelect1.tr")
wgExpSamples <- loadSamples(wgeFN)
wgExpMeans <- rowMeans(as.data.frame(wgExpSamples))
head(wgExpMeans)
}
\keyword{gene expression}
