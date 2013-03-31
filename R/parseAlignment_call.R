

parseAlignment <- function( alignFile, outFile, trSeqFile, inputFormat=NULL, trInfoFile=NULL, expressionFile=NULL, readsN=NULL, uniform=TRUE, limitA=NULL, lenMu=NULL, lenSigma=NULL, verbose=NULL, veryVerbose=NULL, procN=NULL, pretend=FALSE){
   args <- c('parseAlignment', alignFile, '--outFile', outFile , '--trSeqFile', trSeqFile)
   if (!is.null(inputFormat)) {
      args <- c(args, '--format', inputFormat)
   }
   if (!is.null(trInfoFile)) {
      args <- c(args, '--trInfoFile', trInfoFile)
   }
   if (!is.null(expressionFile)) {
      args <- c(args, '--expressionFile', expressionFile)
   }
   if (!is.null(readsN)) {
      args <- c(args, '--readsN', readsN)
   }
   if (!is.null(limitA)) {
      args <- c(args, '--limitA', limitA)
   }
   if (!is.null(lenMu)) {
      args <- c(args, '--lenMu', lenMu)
   }
   if (!is.null(lenSigma)) {
      args <- c(args, '--lenSigma', lenSigma)
   }
   if ((!is.null(uniform)) && (uniform)) {
      args <- c(args, '--uniform')
   }
   if (!is.null(verbose) && (verbose)) {
      args <- c(args, '--verbose')
   }
   if (!is.null(veryVerbose) && (veryVerbose)) {
      args <- c(args, '--veryVerbose')
   }
   if (!is.null(procN)) {
      args <- c(args, '--procN', procN)
   }

   if(pretend){
      writeLines(.specialPaste(args))
   }else{
      argc <- length(args);
      ##dyn.load(paste("src/parseAlignment", .Platform$dynlib.ext, sep=""));
      result <- .C("_parseAlignment", as.integer(argc), as.character(args));
   }
}
