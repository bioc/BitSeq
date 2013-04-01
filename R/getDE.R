## compute PPLR for two or more conditions
getDE <- function( conditions, outPrefix=NULL, samples=FALSE, trInfoFile=NULL, norm=NULL, seed=NULL, pretend=FALSE ){
   ## pretend needs prefix
   if(pretend && is.null(outPrefix)){
      message("In case of no outPrefix provided this function uses R temporary directories which are only valid during the current session.")
      stop("Please provide outPrefix when using the pretend option.")
   }
   ## check file extensions
   ext <- .getExt(conditions[[1]][1]);
   for(it in unlist(conditions)){
      if(!file.exists(it))stop(paste("File",it,"does not exist"));
      if(ext!=.getExt(it))stop(paste("File",it,"has different extension from first file."));
   }
   if(is.null(outPrefix)){
      outPrefix <- tempfile("dataBS-DE-");
   }
   if(is.null(trInfoFile)){
      trInfoFile=.changeExt(conditions[[1]][1]);
      if(file.exists(trInfoFile)){
         message(paste("## No trInfoFile provided, will try using",trInfoFile,"for result's row names."));
      }else{
         trInfoFile=NULL;
      }
   }
   meanFile <- paste(outPrefix, "Lmean", sep=".");
   parFile <- paste(outPrefix, "par", sep=".");
   message("Computing overall mean.");
   getMeanVariance(conditions,meanFile,log=TRUE,norm=norm,pretend=pretend);
   message("Estimating hyperparameters.");
   estimateHyperPar(parFile,conditions,meanFile=meanFile,norm=norm,seed=seed,pretend=pretend);
   message("Estimating condition mean expression and PPLR.");
   estimateDE(conditions,outPrefix,parFile,samples=samples,norm=norm,seed=seed,pretend=pretend);
   if(samples){
      samplesFiles = c();
      for(i in 1:length(conditions)){
         samplesFiles = c(samplesFiles, paste(outPrefix,"-C",(i-1),".est",sep=""));
      }
   }else{
      samplesFiles = NULL;
   }
   pplrFN <- paste(outPrefix, "pplr", sep=".");
   if(pretend){
      data <- NULL
   }else{
      data <- loadSamples( pplrFN, trInfoFile);
      condN <- length(conditions);
      pplrList <- c();
      fcList <- c();
      meanList <- c();
      for(i in 1:condN){
         meanList <- c(meanList, paste("mean",i,sep=" "));
         if(i == condN) break;
         for(j in (i+1):condN){
            pplrList <- c(pplrList, paste("pplr ",i,"~",j,sep=""));
            fcList <- c(fcList, paste("log2FC ",i,"~",j,sep=""), "ciLow", "ciHigh");
         }
      }
      colnames(data)<-c(pplrList, fcList, meanList);
   }
   return(list(pplr=data,fn=list(pplr=pplrFN,samplesFiles=samplesFiles)));
}

