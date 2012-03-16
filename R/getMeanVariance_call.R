#   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
#   args.addOptionB("l","log","log",0,"Use logged values.");
#   args.addOptionS("t","type","type",0,"Type of variance, possible values: [sample,sqDif] for sample variance or squared difference.","sample");

getMeanVariance <- function(sampleFiles, outFile, log=NULL, type=NULL, verbose=NULL){

   args <- c('getMeanVariance', sampleFiles, '--outFile', outFile);
   if ((!is.null(log)) && (log)) {
      args <- c(args, '--log')
   }
   if (!is.null(type)) {
      args <- c(args, '--type', type)
   } 
   if (!is.null(verbose) && (verbose)) {
      args <- c(args, '--verbose')
   }

   argc <- length(args);
   ## dyn.load(paste("src/getVariance", .Platform$dynlib.ext, sep=""));
   result <- .C("_getVariance", as.integer(argc), as.character(args));
}
