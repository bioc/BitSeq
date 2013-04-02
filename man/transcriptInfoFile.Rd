\name{transcriptInfoFile}
\alias{tri.file.hasGeneNames}
\alias{tri.file.setGeneNames}
\title{Manage file containing information about transcript reference}
\description{Manage file containing information about the transcript reference. These functions are used for verifying and updating transcript information DataFrame.}
\usage{
tri.file.hasGeneNames(trInfoFile)
tri.file.setGeneNames(trInfoFile, geneNames, transcriptNames=NULL)
}
\arguments{
\item{trInfoFile}{Name of the file containing transcript information or the where the information should be stored.}
\item{geneNames}{Vector with new gene names that should be assigned to transcripts.}
\item{transcriptNames}{Names of transcripts that should be associated with the gene names.}
}
\details{
If not provided with the information, \code{BitSeq} extracts information about the transcript reference from the alignment and sequence files.
This information is stored in so called transcript information(trInfo) file, usually having extension \code{.tr}.
This file contains columns with gene names (if available), transcript names, transcript lengths and optionally with adjusted lengths of transcripts. 
Important use of trInfo file is for calculating gene expression or within gene expression, where the file is used for determining which transcripts belong to which genes.
However, for this the gene names have to be properly set in the transcript info, which is not always the case.

Function \code{tri.file.hasGeneNames} determines whether gene names are properly set in the transcript information file and returns \code{TRUE} or \code{FALSE} and a warning message identifying the problem.

Function \code{tri.file.setGeneNames} updates the gene names of a transcript information in file provided by argument \code{trInfoFile}. 
The vector \code{geneNames} should provide gene names of transcripts and be of the same length as is the number of transcripts.
The gene names have to be either ordered as their appropriate transcripts in \code{trInfoFile} file, or if ordered differently, vector of transcript names, ordered as gene names has to be provided by argument \code{transcriptNames}. The names in \code{transcriptNames} have to correspond to the transcript names in actual file.
}
\value{
Function \code{tri.file.hasGeneNames} returns boolean value.
}
\author{Peter Glaus}

\seealso{\code{\link{getExpression}}, \code{\link{getGeneExpression}}, \code{\link{tri.load}}, \code{\link{tri.save}}}
\examples{
setwd(system.file("extdata",package="BitSeq"))
## this should be true
tri.file.hasGeneNames("ensSelect1.tr")
}
\keyword{transcript information}
