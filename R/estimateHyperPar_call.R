#  args.addOptionB("V","veryVerbose","veryVerbose",0,"More verbose output.");
#  args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
#  args.addOptionS("p","paramsAllFile","paramsAllFileName",0,"Name of the file to which to store all parameter values generated prior to lowess smoothing.");
#  args.addOptionS("","var","varFileName",0,"Name of the joint variance file.");
#  args.addOptionL("g","gourpsNumber","groupsN",0,"Number of groups of transcript of similar size.",100);
#  args.addOptionL("s","samplesNumber","samplesN",0,"Number of samples generated for each group.",SAMPLES_N);
#  args.addOptionD("l","lambda0","lambda0",0,"Precision scaling parameter lambda0.",0.5);
#  args.addOptionD("","exThreshold","exT",0,"Threshold of lowest expression for which the estimation is done.",-5);
#  args.addOptionB("S","smoothOnly","smoothOnly",0,"Input file contains previously sampled hyperparameters which should smoothed only");
#  args.addOptionD("","lowess-f","lowess-f",0,"Parameter F for lowess smoothing specifying amount of smoothing.",0.2);
#  args.addOptionL("","lowess-steps","lowess-steps",0,"Parameter Nsteps for lowess smoothing specifying number of iterations.",5);
#  args.addOptionB("","force","force",0,"Force smoothing",true);

estimateHyperPar <- function( outFile, conditions=NULL, paramsInFile=NULL, meanFile=NULL, force=TRUE, exThreshold=NULL, lambda0=NULL, paramsAllFile=NULL, smoothOnly=NULL, lowess_f=NULL, lowess_steps=NULL, verbose=NULL, veryVerbose=NULL, norm=NULL, seed=NULL, pretend=FALSE ){
   
   ## unlist norm
   norm <- unlist(norm);
   if (is.null(paramsInFile)){
      args <- c('estimateHyperPar', unlist(conditions[[1]]));
      ## parse other conditions
      for(i in 2:length(conditions)){
         args <- c(args, 'C', unlist(conditions[[i]]));
      }
      args <- c(args , '--outFile', outFile)
   }else{
      if(is.null(smoothOnly) || (!smoothOnly))stop("Please use smoothOnly option if you only want to smooth previously estimated hyperparameters.");
      args <- c('estimateHyperPar', paramsInFile, '--outFile', outFile)
   }
   if (!is.null(meanFile)) {
      args <- c(args, '--meanFile', meanFile)
   }
   if (!is.null(exThreshold)) {
      args <- c(args, '--exThreshold', exThreshold)
   }
   if (!is.null(lambda0)) {
      args <- c(args, '--lambda0', lambda0)
   }
   if (!is.null(paramsAllFile)) {
      args <- c(args, '--paramsAllFile', paramsAllFile)
   }
   if (!is.null(lowess_f)) {
      args <- c(args, '--lowess-f', lowess_f)
   }
   if (!is.null(lowess_steps)) {
      args <- c(args, '--lowess-steps', lowess_steps)
   }
   if ((!is.null(smoothOnly)) && (smoothOnly)) {
      args <- c(args, '--smoothOnly')
   }
   if (force) {
      args <- c(args, '--noforce')
   }
   if (!is.null(verbose) && (verbose)) {
      args <- c(args, '--verbose')
   }
   if (!is.null(veryVerbose) && (veryVerbose)) {
      args <- c(args, '--veryVerbose')
   }
   if (!is.null(norm)) {
      if(length(unlist(conditions)) != length(norm)){
         stop("The number of normalization constants has to match the number of sample files.");
      }
      args <- c(args, '--norm', paste(norm, collapse=","));
   }
   if (!is.null(seed)) {
      args <- c(args, '--seed', seed)
   }

   if(pretend){
      writeLines(.specialPaste(args))
   }else{
      argc <- length(args);
      ## dyn.load(paste("src/estimateHyperPar", .Platform$dynlib.ext, sep=""));
      result <- .C("_estimateHyperPar", as.integer(argc), as.character(args));
   }
}
