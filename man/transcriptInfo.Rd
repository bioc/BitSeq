\name{transcriptInfo}
\alias{tri.load}
\alias{tri.save}
\alias{tri.hasGeneNames}
\alias{tri.setGeneNames}
\title{Manage information about transcript reference}
\description{Manage information about the transcript reference. These functions are used for reading, saving and updating transcript information DataFrame.}
\usage{
tri.load(trInfoFile)
tri.save(trInfo, trInfoFile)
tri.hasGeneNames(trInfo)
tri.setGeneNames(trInfo, geneNames, transcriptNames=NULL)
}
\arguments{
\item{trInfoFile}{Name of the file containing transcript information or the where the information should be stored.}
\item{trInfo}{\code{DataFrame} containing the transcript information.}
\item{geneNames}{Vector with new gene names that should be assigned to transcripts.}
\item{transcriptNames}{Names of transcripts that should be associated with the gene names.}
}
\details{
If not provided with the information, \code{BitSeq} extracts information about the transcript reference from the alignment and sequence files.
This information is stored in so called transcript information(trInfo) file, usually having extension \code{.tr}.
This file contains columns with gene names (if available), transcript names, transcript lengths and optionally with adjusted lengths of transcripts. 
The expression of transcripts is reported in the same order as are the transcripts ordered in the trInfo file, hence it serves as identification of final results.

Other important use of trInfo file is for calculating gene expression or within gene expression, where the file is used for determining which transcripts belong to which genes.
However, for this the gene names have to be properly set in the transcript info, which is not always the case.

Function \code{tri.load} loads transcript information from a file provided by argument \code{trInfoFile} into a \code{DataFrame}. 

Function \code{tri.save} saves transcript information from a \code{DataFrame} provided by \code{trInfo} argument into a file name provided by argument \code{trInfoFile}.

Function \code{tri.hasGeneNames} determines whether gene names are properly set in the transcript information and returns \code{TRUE} or \code{FALSE} and a warning message identifying the problem.

Function \code{tri.setGeneNames} changes gene names of a transcript information \code{trInfo} and retruns new \code{DataFrame} with updated values.
The vector \code{geneNames} should provide gene names of transcripts and be of the same length as is the number of transcripts.
The gene names have to be either ordered as their appropriate transcripts in \code{trInfo} object, or if ordered differently, vector of transcript names, ordered as gene names has to be provided by argument \code{transcriptNames}. The names in \code{transcriptNames} have to correspond to the transcript names in \code{trInfo} object.
}
\value{
Function \code{tri.load} returns \code{DataFrame} with transcript information.

Function \code{tri.hasGeneNames} returns boolean value.

Function \code{tri.setGeneNames} returns \code{DataFrame} with transcript information containing updated gene names (Note: the transcript names do not change.).
}
\author{Peter Glaus}

\seealso{\code{\link{getExpression}}, \code{\link{getGeneExpression}}, \code{\link{tri.file.setGeneNames}}}
\examples{
setwd(system.file("extdata",package="BitSeq"))
trinfo <- tri.load("ensSelect1.tr")
trinfo[1:10,]
## this should be true
tri.hasGeneNames(trinfo)
## reverse the gene order - this will make the information INCORRECT
rev.trinfo <- tri.setGeneNames(trinfo, rev(trinfo[,1]))
rev.trinfo[1:10,]
tri.save(rev.trinfo, "reversed-ensSelect1.tr")
}
\keyword{transcript information}
