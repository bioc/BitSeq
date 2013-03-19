## load transcript names (and information) into 
loadTranscriptNames <- function(fileName){
   if(is.null(fileName) || (!file.exists(fileName)))stop("Please provide valid file name.");
   trNames <- read.table(fileName, sep=" ");
   ## omit last column if its NA
   if(is.na(trNames[1,dim(trNames)[2]])){
      return( IRanges::DataFrame(trNames[,1:dim(trNames)[2]-1]));
   }else{
      return( IRanges::DataFrame(trNames));
   }
}

## check whether transcript info file has gene names set
hasGeneNames <- function(filename){
   trInfo <- loadTranscriptNames(filename);
   if(is.null(trInfo))return(FALSE);
   geneN <- unique(trInfo[,1]);
   if(length(geneN) == 1){
      warning("There seems to be just one gene name.");
      return(FALSE);
   }
   if(length(geneN) == length(trInfo[,1])){
      warning("There seems to be just one transcript for every gene.");
      return(FALSE);
   }
   return(TRUE);
}

## set gene names in transcript info file
setGeneNames <- function(filename, geneNames, transcriptNames=NULL){
   trInfo <- loadTranscriptNames(filename);
   if(dim(trInfo)[1] != length(geneNames)){
      stop("Number of gene names does not match number of transcripts.\n Please provide one gene name per transcript");
   }
   if(is.null(transcriptNames)){
      for(i in 1:length(geneNames)){
         print(trInfo[i,]);
         trInfo[i,1] <- geneNames[i];
         print(trInfo[i,]);
      }
      saveTranscriptInfo(trInfo, filename);
   }else{
      if(length(geneNames) != length(transcriptNames)){
         stop("Number of gene names does not match number of transcript names.\n Please provide one gene name per transcript name.");
      }
      errorN <- 0;
      for(i in 1:dim(trInfo)[1]){
         ind <- which(transcriptNames == trInfo[i,2]);
         if(length(ind)==0){
            warning(paste("Couldn't find gene name for transcript",trInfo[i,2]));
            errorN <- errorN+1;
         }else if(length(ind)>1){
            warning(paste("Multiple gene names for transcript",trInfo[i,2]));
            errorN <- errorN+1;
         }else {
            message(sprintf("INDEX %d\n",ind));
            trInfo[i,1] <- geneNames[ind];
         }
         if(errorN>4)stop("Too many errors.");
      }
      saveTranscriptInfo(trInfo, filename);
   }
}

## Save transcript information data from DataFrame into file.
saveTranscriptInfo <- function(trInfo, filename){
   print("====");
   stopifnot( is(trInfo,"DataFrame"));
   if(file.exists(filename)){
      file.remove(filename);
   }
   header <- sprintf("# M %i",dim(trInfo)[1]);
   outF <- file(filename, "w");
   writeLines(header,outF);
   message("lines");
   for(i in 1:dim(trInfo)[1]){
      if(length(trInfo[i,])==3){
   message("line");
         line <- sprintf("%s %s %i",trInfo[i,1], trInfo[i,2], trInfo[i,3]);
      }else{
         line <- sprintf("%s %s %i %f",trInfo[i,1], trInfo[i,2], trInfo[i,3], trInfo[i,4]);
      }
      writeLines(line, outF);
   }
   close(outF);
   ##write.table(as.data.frame(trInfo),filename,append=TRUE,row.names=FALSE,col.names=FALSE);
}
