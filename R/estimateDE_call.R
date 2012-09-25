#  args.addOptionS("o","outFile","outFilePrefix",1,"Prefix for the output files.");
#  args.addOptionS("p","parameters","parFileName",1,"File containing estimated hyperparameters.");
#  args.addOptionB("s","samples","samples",0,"Produce samples of condition mean expression apart from PPLR and confidence.");
#  args.addOptionD("l","lambda0","lambda0",0,"Parameter lambda_0.",LAMBDA_0);
#  args.addOptionD("c","confidencePerc","cf",0,"Percentage for confidence intervals.", 5);

estimateDE <- function( cond1, cond2, outFile, parFile, lambda0=NULL, samples=NULL, confidencePerc=NULL, verbose=NULL, pretend=FALSE ){
   
   args <- c('estimateDE',cond1,'C',cond2, '--outPrefix', outFile , '--parameters', parFile)
   if (!is.null(lambda0)) {
      args <- c(args, '--lambda0', lambda0)
   }
   if (!is.null(samples) && (samples)) {
      args <- c(args, '--samples')
   }
   if (!is.null(confidencePerc)) {
      args <- c(args, '--confidencePerc', confidencePerc)
   }
   
   if (!is.null(verbose) && (verbose)) {
      args <- c(args, '--verbose')
   }

   if(pretend){
      writeLines(.specialPaste(args))
   }else{
      argc <- length(args);
      ##dyn.load(paste("src/estimateDE", .Platform$dynlib.ext, sep=""))
      result <- .C("_estimateDE", as.integer(argc), as.character(args))
   }
}
