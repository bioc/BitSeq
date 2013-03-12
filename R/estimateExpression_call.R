## see list of possible command line options at the end.

## just for checking validity of arguments passed to eE from getExpression
.argCheck.eELegacy <- function(probFile, outFile, parFile=NULL, outputType=NULL, gibbs=NULL, trInfoFile=NULL, thetaActFile=NULL, MCMC_burnIn=NULL, MCMC_samplesN=NULL, MCMC_samplesSave=NULL, MCMC_samplesNmax=NULL, MCMC_chainsN=NULL, MCMC_scaleReduction=NULL, MCMC_dirAlpha=NULL, seed=NULL, verbose=NULL, pretend=FALSE){
   return(TRUE);
}
.argCheck.eE <- function(probFile, outFile, parFile=NULL, outputType=NULL, gibbs=NULL, trInfoFile=NULL, thetaActFile=NULL, MCMC_burnIn=NULL, MCMC_samplesN=NULL, MCMC_samplesSave=NULL, MCMC_chainsN=NULL, MCMC_dirAlpha=NULL, seed=NULL, verbose=NULL, procN=NULL, pretend=FALSE){
   return(TRUE);
}


## Call of estimateExpression which uses [new] default method for convergence checking.
estimateExpression <- function(probFile, outFile, parFile=NULL, outputType=NULL, gibbs=NULL, trInfoFile=NULL, thetaActFile=NULL, MCMC_burnIn=NULL, MCMC_samplesN=NULL, MCMC_samplesSave=NULL, MCMC_chainsN=NULL, MCMC_dirAlpha=NULL, seed=NULL, verbose=NULL, procN=NULL, pretend=FALSE){
   args <- c('estimateExpression ',probFile, '--outPrefix', outFile)
   if (!is.null(parFile)){
      args <- c(args, '--parFile', parFile )
   }
   if (!is.null(outputType)) {
      args <- c(args, '--outType', outputType)
   }
   if ((!is.null(gibbs)) && (gibbs)){
        args <- c(args, '--gibbs')
   }
   if (!is.null(trInfoFile)) {
      args <- c(args, '--trInfoFile', trInfoFile )
   }
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
   if (!is.null( MCMC_chainsN)) {
      args <- c(args, '--MCMC_chainsN', MCMC_chainsN)
   }
   if (!is.null( MCMC_dirAlpha)) {
      args <- c(args, '--MCMC_dirAlpha', MCMC_dirAlpha)
   }
   if (!is.null(seed)) {
      args <- c(args, '--seed', seed)
   }
   if (!is.null(verbose) && (verbose)) {
      args <- c(args, '--verbose')
   }
   if (!is.null(procN)) {
      args <- c(args, '--procN', procN)
   }

   ## print(args)

   if(pretend){
      writeLines(.specialPaste(args))
   }else{
      argc <- length(args);
      ## dyn.load(paste("src/estimateExpression", .Platform$dynlib.ext, sep=""));
      result <- .C('_estimateExpression', as.integer(argc), as.character(args));
   }
}

## Legacy call of estimateExpression which uses old method for convergence checking based on "--scaleReduction".
estimateExpressionLegacy <- function(probFile, outFile, parFile=NULL, outputType=NULL, gibbs=NULL, trInfoFile=NULL, thetaActFile=NULL, MCMC_burnIn=NULL, MCMC_samplesN=NULL, MCMC_samplesSave=NULL, MCMC_samplesNmax=NULL, MCMC_chainsN=NULL, MCMC_scaleReduction=NULL, MCMC_dirAlpha=NULL, seed=NULL, verbose=NULL, pretend=FALSE){
## , procN=NULL
   args <- c('estimateExpression', '--scaleReduction',probFile, '--outPrefix', outFile)
   if (!is.null(parFile)){
      args <- c(args, '--parFile', parFile )
   }
   if (!is.null(outputType)) {
      args <- c(args, '--outType', outputType)
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
   if (!is.null(seed)) {
      args <- c(args, '--seed', seed)
   }
   if (!is.null(verbose) && (verbose)) {
      args <- c(args, '--verbose')
   }
   ## print(args)

   if(pretend){
      writeLines(.specialPaste(args))
   }else{
      argc <- length(args);
      ## dyn.load(paste("src/estimateExpression", .Platform$dynlib.ext, sep=""));
      result <- .C('_estimateExpression', as.integer(argc), as.character(args));
   }
}

## Options:
##   --help
##     Show this help information.
## 
##   --MCMC_burnIn=<MCMC_burnIn>
##     Length of sampler's burn in period. (default: 1000)
## 
##   --MCMC_chainsN=<MCMC_chainsN>
##     Number of parallel chains used. At least two chains will be used. (default: 4)
## 
##   --MCMC_dirAlpha=<MCMC_dirAlpha>
##     Alpha parameter for the Dirichlet distribution. (default: 1)
## 
##   --MCMC_samplesDOmax
##     Produce maximum number of samples (samplesNmax) in second iteration and quit. (default: Off)
## 
##   --MCMC_samplesN=<MCMC_samplesN>
##     Initial number of samples produced. Doubles after every iteration. (default: 1000)
## 
##   --MCMC_samplesNmax=<MCMC_samplesNmax>
##     Maximum number of samples produced in one iteration. After producing samplesNmax samples sampler finishes. (default: 50000)
## 
##   --MCMC_samplesSave=<MCMC_samplesSave>
##     Number of samples recorder in total. (default: 1000)
## 
##   --MCMC_scaleReduction=<MCMC_scaleReduction>
##     Target scale reduction, sampler finishes after this value is met. (default: 1.2)
## 
##   -G ,   --gibbs
##     Use gibbs sampling instead of collapsed gibbs sampling. (default: Off)
## 
##   -o <outFilePrefix> ,   --outPrefix=<outFilePrefix>
##     Prefix for the output files.
## 
##   -O <outputType> ,   --outType=<outputType>
##     Output type (theta, RPKM, counts, tau). (default: theta)
## 
##   -p <parFileName> ,   --parFile=<parFileName>
##     File containing parameters for the sampler, which can be otherwise specified by --MCMC* options. As the file is checked after every MCMC iteration, the parameters can be adjusted while running.
## 
##   -P <procN> ,   --procN=<procN>
##     Limit the maximum number of threads to be used. (Default is the number of MCMC chains.)
## 
##   --scaleReduction
##     Use scale reduction as stopping criterion, instead of computing effective sample size. (default: Off)
## 
##   -s <seed> ,   --seed=<seed>
##     Random initialization seed.
## 
##   --thetaActFile=<thetaActFileName>
##     File for logging noise parameter theta^{act}.
## 
##   -t <trInfoFileName> ,   --trInfoFile=<trInfoFileName>
##     File containing transcript information. (Necessary for RPKM)
## 
##   -v ,   --verbose
##     Verbose output. (default: Off)
