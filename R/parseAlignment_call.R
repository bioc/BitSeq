#   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
#   args.addOptionS("f","format","format",0,"Input format: either SAM, BAM or bowtie's MAP[not implemented yet].","SAM");
#   args.addOptionS("t","trInfoFile","trInfoFileName",0,"If transcript(reference sequence) information is contained within SAM file, program will write this information into <trInfoFile>, otherwise it will look for this information in <trInfoFile>.");
#   args.addOptionS("s","trSeqFile","trSeqFileName",1,"Transcript sequence in FASTA format --- for non-uniform read distribution estimation.");
#   args.addOptionS("e","expressionFile","expFileName",0,"Transcript relative expression estimates --- for better non-uniform read distribution estimation.");
#   args.addOptionL("N","readsN","readsN",0,"Total number of reads. This is not necessary if [SB]AM contains also reads with no valid alignments.");
#   args.addOptionS("","failed","failed",0,"File name where to save names of reads that failed to align as pair.");
#   args.addOptionB("","uniform","uniform",0,"Use uniform read distribution.");
#   args.addOptionD("","lenMu","lenMu",0,"Set mean of log fragment length distribution. (l_frag ~ LogNormal(mu,sigma^2)");
#   args.addOptionD("","lenSigma","lenSigma",0,"Set sigma^2 (or variance) of log fragment length distribution. (l_frag ~ LogNormal(mu,sigma^2)");
#   args.addOptionS("","distributionFile","distributionFileName",0,"Name of file to which read-distribution should be saved.");
#   args.addOptionL("P","procN","procN",0,"Maximum number of threads to be used.",8);
#   args.addOptionB("V","veryVerbose","veryVerbose",0,"Very verbose output.");
#   args.addOptionL("","noiseMismatches","numNoiseMismatches",0,"Number of mismatches to be considered as noise.",LOW_PROB_MISSES);


parseAlignment <- function( alignFile, outFile, trSeqFile, inputFormat=NULL, trInfoFile=NULL, expressionFile=NULL, readsN=NULL, uniform=TRUE, lenMu=NULL, lenSigma=NULL, verbose=NULL, veryVerbose=NULL, pretend=FALSE){
##  procN=NULL,
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
   ##if (!is.null(procN)) {
   ##   args <- c(args, '--procN', procN)
   ##}
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

   if(pretend){
      print(paste(args,collapse=" "))
   }else{
      argc <- length(args);
      ##dyn.load(paste("src/parseAlignment", .Platform$dynlib.ext, sep=""));
      result <- .C("_parseAlignment", as.integer(argc), as.character(args));
   }
}
