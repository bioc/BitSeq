\name{getGeneExpression}
\alias{getGeneExpression}
\alias{getWithinGeneExpression}
\title{Calculate gene expression or relative within gene expression}
\description{Calculate either gene expression or relative within gene expression using transcript expression samples and transcript information file.}
\usage{
getGeneExpression(sampleFile, outFile=NULL, trInfoFile=NULL)
getWithinGeneExpression(sampleFile, outFile=NULL, trInfoFile=NULL)
}
\arguments{
\item{sampleFile}{File containing the transcript expression samples.}
\item{outFile}{Name of the output file. If not used, function uses temporary file.}
\item{trInfoFile}{Transcript information file. If not used, function tries file with same name and extension \code{tr}. The file has to contain valid gene transcript mapping, see detail below.}
}
\details{
The \code{getGeneExpression} function takes samples of transcript expression and produces file with expression of genes by adding up transcript expression.

The \code{getWithinGeneExpression} function takes samples of transcript expression and produces file with relative within gene expression samples for each transcript.

Both function need valid transcript information file which contains gene transcript mapping. 
The first line should contain "# M <numberOfTranscripts>" and the following \code{numberOfTranscripts} lines have to contain "<geneName> <transcriptName> <transcriptLength>". 
Example is provided in \code{extdata/ensSelect1.tr}. 
Please note that the transcript information file automatically generated from alignment files are not sufficient because SAM/BAM files do not include gene names. 
We hope to provide more convenient way in future versions of BitSeq.
}
\value{
Name of file containing the new expression samples.
}
\author{Peter Glaus}

\seealso{\code{\link{getExpression}}}
\examples{
setwd(system.file("extdata",package="BitSeq"))
## gene expression
getGeneExpression("data-c0b1.rpkm", "data-c0b1-GE.rpkm", "ensSelect1.tr")
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
