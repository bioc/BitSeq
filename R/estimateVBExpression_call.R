
estimateVBExpression <- function(probFile, outFile, outputType=NULL, trInfoFile=NULL, seed=NULL, samples=NULL, optLimit=1e-5, optMethod="FR", procN=4, verbose=FALSE, veryVerbose=FALSE, pretend=FALSE){
   args <- c('estimateExpression ',probFile, '--outPrefix', outFile)
   if (!is.null(outputType)) {
      args <- c(args, '--outType', outputType)
   }
   if (!is.null(trInfoFile)) {
      args <- c(args, '--trInfoFile', trInfoFile )
   }
   if (!is.null(seed)) {
      args <- c(args, '--seed', seed)
   }
   if (!is.null( samples)) {
      args <- c(args, '--samples',samples )
   }
   if (!is.null( optLimit)) {
      args <- c(args, '--optLimit',optLimit )
   }
   if (!is.null( optMethod)) {
      args <- c(args, '--method',optMethod )
   }
   if (!is.null(procN)) {
      args <- c(args, '--procN', procN)
   }
   if (!is.null(verbose) && (verbose)) {
      args <- c(args, '--verbose')
   }
   if (!is.null(veryVerbose) && (veryVerbose)) {
      args <- c(args, '--veryVerbose')
   }

   ## print(args)

   if(pretend){
      writeLines(.specialPaste(args))
   }else{
      argc <- length(args);
      result <- .C('_estimateVBExpression', as.integer(argc), as.character(args));
   }
}

