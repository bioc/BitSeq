## compute transcript expression
getExpression <- function(alignFile, trSeqFile, outPrefix=NULL, uniform=TRUE, type="RPKM", log=FALSE, limitA=NULL, seed=NULL, pretend=FALSE, ... ){
   if(pretend && is.null(outPrefix)){
      message("In case of no outPrefix provided this function uses R temporary directories which are only valid during the current session.")
      stop("Please provide outPrefix when using the pretend option.")
   }
   if(is.null(alignFile) || (!file.exists(alignFile)))stop("Please provide valid sam/bam file as alignFile.");
   if(is.null(trSeqFile) || (!file.exists(trSeqFile)))stop("Please provide valid Fasta reference file as trSeqFile.");
   ## check type:
   type <- tolower(type);
   types <- c("rpkm","theta","counts");
   if(sum(types == type) != 1)stop(paste("Unknown output type: ", type));
   ## set prefix
   if(is.null(outPrefix)){
      fName <- tail(unlist(strsplit(tail( unlist(strsplit(.removeExt(alignFile),"/")),1),"\\\\")),1)
      if(fName != ""){
         outPrefix <- tempfile(paste(fName,"-BS-",sep=""));
      }else{
         outPrefix <- tempfile("dataBS-E-");
      }
   }
   probF <- paste(outPrefix,"prob",sep=".");
   trF <- paste(outPrefix,"tr",sep=".");
   ## check arguments that will be passed to estimateExpression later:
   .argCheck.eE(probF, outPrefix, outputType=type, trInfoFile=trF, ... );

   if(! uniform){
      probUF <- paste(outPrefix,"-U",".prob",sep="");
      trUF <- paste(outPrefix,"-U",".tr",sep="");
      message("## Pre-Computing alignment probabilities with uniform read distribution.");
      parseAlignment(alignFile, probUF, trSeqFile, trInfoFile=trUF, uniform=TRUE, limitA=limitA, pretend=pretend);
      message("## Pre-Computing expression with uniform read distribution.");
      if (!is.null(seed)) {
         # Use different seed for initial sampling.
         seed <- seed + 17;
      }
      estimateVBExpression(probUF, paste(outPrefix,"-U",sep=""), trInfoFile=trUF, seed=seed, pretend=pretend);
      exprFile <- paste(outPrefix,"-U",".m_alphas",sep=""); 
   }else{
      exprFile <- NULL;
   }
   
   message("## Computing alignment probabilities.");
   parseAlignment(alignFile, probF, trSeqFile, trInfoFile=trF, uniform=uniform, limitA=limitA, expressionFile=exprFile,pretend=pretend);
   message("## Estimating transcript expression levels.");
   estimateExpression(probF, outPrefix, outputType=type, trInfoFile=trF, seed=seed, pretend=pretend, ... );
   ## Generated files:
   outFile <- paste(outPrefix, type, sep=".");
   thetaMeansFile <- paste(outPrefix, "thetaMeans", sep=".");
   meanFile <- paste(outPrefix, "mean", sep=".");
   message("## Computing means.");
   getMeanVariance(c(outFile), meanFile, log=log, pretend=pretend );
   ## load means, trInfo and counts;
   if(pretend){
      means <- NULL;
      trInfo <- NULL;
      counts <- NULL;
   }else{
      means <- loadSamples( meanFile )
      means[,2] <- sqrt(means[,2])
      colnames(means)<-c("mean","stdev");
      counts <- read.table(thetaMeansFile, sep=" ")[,3]
      trInfo <- tri.load(trF);
   }
   return(list(exp=means , fn=outFile, counts=counts, trInfo=trInfo));
}

