## compute PPLR for two conditions
getDE <- function( cond1, cond2, outPrefix=NULL, samples=FALSE, trInfoFile=NULL, pretend=FALSE ){
   ## pretend needs prefix
   if(pretend && is.null(outPrefix)){
      message("In case of no outPrefix provided this function uses R temporary directories which are only valid during the current session.")
      stop("Please provide outPrefix when using the pretend option.")
   }
   ## check file extensions
   ext <- .getExt(cond1[1]);
   for(it in c(cond1,cond2)){
      if(!file.exists(it))stop(paste("File",it,"does not exist"));
      if(ext!=.getExt(it))stop(paste("File",it,"has different extension from first file."));
   }
   if(is.null(outPrefix)){
      outPrefix <- tempfile("dataBS-DE-");
   }
   if(is.null(trInfoFile)){
      trInfoFile=.changeExt(cond1[1]);
      if(file.exists(trInfoFile)){
         message(paste("## No trInfoFile provided, will try using",trInfoFile,"for result's row names."));
      }else{
         trInfoFile=NULL;
      }
   }
   meanFile <- paste(outPrefix, "Lmean", sep=".");
   parFile <- paste(outPrefix, "par", sep=".");
   message("Computing overall mean.");
   getMeanVariance(c(cond1,cond2),meanFile,log=TRUE,pretend=pretend);
   message("Estimating hyperparameters.");
   estimateHyperPar(parFile,cond1,cond2,meanFile=meanFile,pretend=pretend);
   message("Estimating condition mean expression and PPLR.");
   estimateDE(cond1,cond2,outPrefix,parFile,samples=samples,pretend=pretend);
   if(samples){
      c1Res <- paste(outPrefix,"-C0.est",sep="");
      c2Res <- paste(outPrefix,"-C1.est",sep="");
   }else{
      c1Res<-NULL;
      c2Res<-NULL;
   }
   pplrFN <- paste(outPrefix, "pplr", sep=".");
   if(pretend){
      data <- NULL
   }else{
      data <- loadSamples( pplrFN, trInfoFile);
      colnames(data)<-c("pplr", "ConfidenceLow", "ConfidenceHigh", "log2FC", "meanC1", "meanC2");
   }
   return(list(pplr=data,fn=list(pplr=pplrFN,C1samples=c1Res,C2samples=c2Res)));
}

## compute transcript expression
getExpression <- function(alignFile, trSeqFile, outPrefix=NULL, uniform=TRUE, type="RPKM", log=FALSE, pretend=FALSE, ... ){
   if(pretend && is.null(outPrefix)){
      message("In case of no outPrefix provided this function uses R temporary directories which are only valid during the current session.")
      stop("Please provide outPrefix when using the pretend option.")
   }
   if(is.null(alignFile) || (!file.exists(alignFile)))stop("Please provide valid sam/bam file as alignFile.");
   if(is.null(trSeqFile) || (!file.exists(trSeqFile)))stop("Please provide valid Fasta reference file as trSeqFile.");
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

   iFormat <- switch(.getExt(alignFile),bam="BAM",BAM="BAM","SAM");
   if(! uniform){
      probUF <- paste(outPrefix,"-U",".prob",sep="");
      print(probUF);
      message("## Pre-Computing alignment probabilities with uniform read distribution.");
      parseAlignment(alignFile, probUF, trSeqFile, inputFormat=iFormat, uniform=FALSE,pretend=pretend);
      message("## Pre-Computing expression with uniform read distribution.");
      estimateExpression(probUF, paste(outPrefix,"-U",sep=""), outputType="theta", MCMC_burnIn=1000, MCMC_samplesN=1000, MCMC_samplesSave=10, MCMC_chainsN=2, pretend=pretend);
      exprFile <- paste(outPrefix,"-U",".thetaMeans",sep=""); 
   }else{
      exprFile <- NULL;
   }
   
   message("## Computing alignment probabilities.");
   parseAlignment(alignFile, probF, trSeqFile, trInfoFile=trF, inputFormat=iFormat, uniform=uniform, expressionFile=exprFile,pretend=pretend);
   message("## Estimating transcript expression levels.");
   estimateExpression(probF, outPrefix, outputType=type, trInfoFile=trF, pretend=pretend, ... );
   message("## Computing means.");
   if(type=="RPKM")type<-"rpkm";
   outFile <- paste(outPrefix, type, sep=".");
   meanFile <- paste(outPrefix, "mean", sep=".");
   getMeanVariance(c(outFile), meanFile, log=log, pretend=pretend );
   if(pretend){
      means <- NULL
   }else{
      means <- loadSamples( meanFile )
      means[,1] <- sqrt(means[,1])
      colnames(means)<-c("mean","stdev");
   }
   return(list(exp=means , fn=outFile));
}


## load expression samples into DataFrame
loadSamples <- function(fileName, trInfoFile=NULL){
   if(is.null(fileName) || (!file.exists(fileName)))stop("Please provide valid file name.");
   trNames <- NULL;
   if(! is.null(trInfoFile)){
      if(file.exists(trInfoFile))
         trNames <- read.table(trInfoFile, sep=" ")[[2]];
   }else{
      trInfoFile <- .changeExt(fileName);
      if(file.exists(trInfoFile)){
         trNames <- read.table(trInfoFile, sep=" ")[[2]];
      }
   }
   data <- read.table(fileName, sep=" ");
   if(is.na(data[1,dim(data)[2]])){
      ## omit last column if its NA
      return( IRanges::DataFrame(data[,1:dim(data)[2]-1], row.names=as.vector(trNames) ));
   }else{
      return( IRanges::DataFrame(data, row.names=as.vector(trNames) ));
   }
}

## write samples into file
writeSamples <- function(data,fileName){
   stopifnot( is(data,"DataFrame"));
   if(file.exists(fileName)){
      file.remove(fileName);
   }
   header <- sprintf("# T (transposed)\n# M %i\n# N %i",dim(data)[1],dim(data)[2]);
   writeLines(header,fileName);
   write.table(as.data.frame(data),fileName,append=TRUE,row.names=FALSE,col.names=FALSE);
}


# change extension of filename
.changeExt <- function(name, ext="tr"){
   nameSplit <- unlist(strsplit(name,"\\."));
   newName <- paste(paste(head(nameSplit,-1), collapse="."), ext, sep=".", collapse=".");
   return(newName);
}
# remove extension of filename
.removeExt <- function(name){
   return(paste(head( unlist(strsplit(name,"\\."))  ,-1), collapse="."));
}
# get extension of filename
.getExt <- function(name){
   return(tail(unlist(strsplit(name,split="\\.")),1));
}

# paste string so that is within 75ch width and lines are ending with \ (i.e. it's usable in bash
.specialPaste <-function(strVec){
   ret <- strVec[1];
   c_len <- nchar(ret);
   for(str in tail(strVec, -1)){
      if(c_len+nchar(str)<75){
         ret <- paste(ret, str);
         c_len <- c_len + nchar(str) + 1;
      }else{
         ret <- sprintf("%s\\\n %s",ret,str);
         c_len <- nchar(str)+1;
      }
   }
   return(ret);
}
