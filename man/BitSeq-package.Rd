\name{BitSeq-package}
\alias{BitSeq-package}
\alias{BitSeq}
\docType{package}
\title{Bayesian Inference of Transcripts from Sequencing data}
\description{
The BitSeq package is targeted for transcript expression analysis and differential expression analysis of RNA-seq data in two stage process. In the first stage it uses Bayesian inference methodology to infer expression of individual transcripts from individual RNA-seq experiments. The second stage of BitSeq embraces the differential expression analysis of transcript expression. Providing expression estimates from replicates of multiple conditions, Log-Normal model of the estimates is used for inferring the condition mean transcript expression and ranking the transcripts based on the likelihood of differential expression.
}
\details{
\tabular{ll}{
Package: \tab BitSeq\cr
Type: \tab Package\cr
Version: \tab 0.3.0\cr
Date: \tab 2012-03-09\cr
License: \tab Artistic-2.0 + other\cr
}
For details of using the package please refer to the Vignette.
}
\author{
Peter Glaus, Antti Honkela and Magnus Rattray
Maintainer: Peter Glaus <glaus@cs.man.ac.uk>
}
\references{
Glaus, P., Honkela, A. and Rattray, M. (2012). Identifying differentially expressed transcripts from RNA-seq data with biological variation. Bioinformatics, 28(13), 1721-1728.
}
\examples{\dontrun{
## basic use
res1 <- getExpression("data-c0b0.sam","ensSelect1.fasta")
res2 <- getExpression("data-c0b1.sam","ensSelect1.fasta")
res3 <- getExpression("data-c1b0.sam","ensSelect1.fasta")
res4 <- getExpression("data-c1b1.sam","ensSelect1.fasta")

deRes <- getDE( list(c(res1$fn, res2$fn), 
                     c(res3$fn, res4$fn)) )
## top 10 differentially expressed
head(deRes$pplr[ order(abs(0.5-deRes$pplr$pplr), decreasing=TRUE ), ], 10)

## advanced use, keeping the intermediate files
parseAlignment( "data-c0b0.sam", 
   outFile = "data-c0b0.prob", 
   trSeqFile = "ensSelect1.fasta",
   trInfoFile = "data.tr",
   uniform = TRUE,
   verbose = TRUE )

estimateExpression( "data-c0b0.prob", 
   outFile = "data-c0b0", 
   outputType = "RPKM",
   trInfoFile = "data.tr",
   MCMC_burnIn = 200,
   MCMC_samplesN = 200,
   MCMC_samplesSave = 100,
   MCMC_scaleReduction = 1.1,
   MCMC_chainsN = 2 )

cond1Files = c("data-c0b0.rpkm","data-c0b1.rpkm")
cond2Files = c("data-c1b1.rpkm","data-c1b1.rpkm")
allConditions = list(cond1Files, cond2Files)

getMeanVariance( allConditions, 
   outFile = "data.means",
   log = TRUE )

estimateHyperPar( allConditions,
   outFile = "data.par",
   meanFile = "data.means",
   verbose = TRUE )

estimateDE( allConditions,
   outFile = "data",
   parFile = "data.par" )

}}
\keyword{ package }
%\seealso{%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~%%}
