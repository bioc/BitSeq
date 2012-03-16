#  args.addOptionS("o","outFile","outFilePrefix",1,"Prefix for the output files.");
#  args.addOptionS("O","outputType","outputType",0,"Output type (theta, RPKM, COV[ERAGE], tau).","theta");
#  args.addOptionB("G","gibbs","gibbs",0,"Use gibbs sampling instead of collapsed gibbs sampling.");
#  args.addOptionS("p","parFile","parFileName",0,"File containing parameters for the sampler, which can be otherwise specified by --MCMC* options. As the file is checked after every MCMC iteration, the parameters can be adjusted while running.");
#  args.addOptionS("t","trInfoFile","trInfoFileName",0,"File containing transcript information. (Necessary for RPKM)");
#  args.addOptionL("P","procN","procN",0,"Maximum number of threads to be used.",8);
#  args.addOptionS("","thetaActFile","thetaActFileName",0,"File for logging noise parameter \theta^{act}.");
#  args.addOptionL("","MCMC_burnIn","MCMC_burnIn",0,"Length of sampler's burn in period.",1000);
#  args.addOptionL("","MCMC_samplesN","MCMC_samplesN",0,"Initial number of samples produced. Doubles after every iteration.",1000);
#  args.addOptionL("","MCMC_samplesSave","MCMC_samplesSave",0,"Number of samples recorder for each chain at the end.",500);
#  args.addOptionL("","MCMC_samplesNmax","MCMC_samplesNmax",0,"Maximum number of samples produced in one iteration. After producing samplesNmax samples sampler finishes.",50000);
#  args.addOptionL("","MCMC_chainsN","MCMC_chainsN",0,"Number of parallel chains used. At least two chains will be used.",4);
#  args.addOptionD("","MCMC_scaleReduction","MCMC_scaleReduction",0,"Target scale reduction, sampler finishes after this value is met.",1.2);
#  args.addOptionD("","MCMC_dirAlpha","MCMC_dirAlpha",0,"Alpha parameter for the Dirichlet distribution.",1.0);


## just for checking validity of arguments passed to eE from getExpression
.argCheck.eE <- function(probFile, outFile, parFile=NULL, outputType=NULL, gibbs=NULL, trInfoFile=NULL, thetaActFile=NULL, MCMC_burnIn=NULL, MCMC_samplesN=NULL, MCMC_samplesSave=NULL, MCMC_samplesNmax=NULL, MCMC_chainsN=NULL, MCMC_scaleReduction=NULL, MCMC_dirAlpha=NULL, verbose=NULL){
   return(TRUE);
}

estimateExpression <- function(probFile, outFile, parFile=NULL, outputType=NULL, gibbs=NULL, trInfoFile=NULL, thetaActFile=NULL, MCMC_burnIn=NULL, MCMC_samplesN=NULL, MCMC_samplesSave=NULL, MCMC_samplesNmax=NULL, MCMC_chainsN=NULL, MCMC_scaleReduction=NULL, MCMC_dirAlpha=NULL, verbose=NULL){
## , procN=NULL
   args <- c('estimateExpression',probFile, '--outFile', outFile)
   if (!is.null(parFile)){
      args <- c(args, '--parFile', parFile )
   }
   if (!is.null(outputType)) {
      args <- c(args, '--outputType', outputType)
   }
   if ((!is.null(gibbs)) && (gibbs)){
        args <- c(args, '--gibbs')
   }
   if (!is.null(trInfoFile)) {
      args <- c(args, '--trInfoFile', trInfoFile )
   }
   ##if (!is.null(procN)) {
   ##   args <- c(args, '--procN', procN)
   ##}
   if (!is.null(thetaActFile)) {
      args <- c(args, '--thetaActFile', thetaActFile)
   }
   if (!is.null( MCMC_burnIn)) {
      args <- c(args, '--MCMC_burnIn',MCMC_burnIn )
   }
   if (!is.null( MCMC_samplesN)) {
      args <- c(args, '--MCMC_samplesN', MCMC_samplesN)
   }
   if (!is.null( MCMC_samplesSave)) {
      args <- c(args, '--MCMC_samplesSave',MCMC_samplesSave )
   }
   if (!is.null( MCMC_samplesNmax)) {
      args <- c(args, '--MCMC_samplesNmax', MCMC_samplesNmax)
   }
   if (!is.null( MCMC_chainsN)) {
      args <- c(args, '--MCMC_chainsN', MCMC_chainsN)
   }
   if (!is.null( MCMC_scaleReduction)) {
      args <- c(args, '--MCMC_scaleReduction', MCMC_scaleReduction)
   }
   if (!is.null( MCMC_dirAlpha)) {
      args <- c(args, '--MCMC_dirAlpha', MCMC_dirAlpha)
   }
   if (!is.null(verbose) && (verbose)) {
      args <- c(args, '--verbose')
   }
   ## print(args)

   argc <- length(args);
   ## dyn.load(paste("src/estimateExpression", .Platform$dynlib.ext, sep=""));
   result <- .C('_estimateExpression', as.integer(argc), as.character(args));
}
