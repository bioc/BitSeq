## load transcript names (and information) into result
tri.load <- function(fileName){
   if(is.null(fileName) || (!file.exists(fileName)))stop("Please provide valid file name.");
   trNames <- read.table(fileName, sep=" ", as.is=c(1,2));

   ## omit last column if its NA
   if(is.na(trNames[1,dim(trNames)[2]])){
      ret <- IRanges::DataFrame(trNames[,1:dim(trNames)[2]-1]);
   }else{
      ret <- IRanges::DataFrame(trNames);
   }
   ## set column names
   if(dim(ret)[2] == 3){
      colnames(ret) <- c("gene", "transcript", "length");
   }else{
      colnames(ret) <- c("gene", "transcript", "length", "adjusted length");
   }
   return(ret);
}

## check whether transcript info has gene names set
tri.hasGeneNames <- function(trInfo){
   if(is.null(trInfo))return(FALSE);
   geneN <- unique(trInfo[,1]);
   if(length(geneN) == 1){
      warning("There seems to be just one gene.");
      return(FALSE);
   }
   if(length(geneN) == length(trInfo[,1])){
      warning("There seems to be just one transcript for every gene.");
      return(FALSE);
   }
   return(TRUE);
}

## check whether transcript info file has gene names set
tri.file.hasGeneNames <- function(filename){
   trInfo <- tri.load(filename);
   return(tri.hasGeneNames(trInfo));
}

## update DataFrame containing transcript info with gene names
tri.setGeneNames <- function(trInfo, geneNames, transcriptNames=NULL){
   stopifnot( is(trInfo,"DataFrame"));
   if(dim(trInfo)[1] != length(geneNames)){
      stop("Number of gene names does not match number of transcripts.\n Please provide one gene name per transcript");
   }
   if(is.null(transcriptNames)){
      for(i in 1:length(geneNames)){
         trInfo[i,1] <- geneNames[i];
      }
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
            trInfo[i,1] <- geneNames[ind];
         }
         if(errorN>4)stop("Too many errors.");
      }
   }
   return(trInfo);
}

## set gene names in transcript info file
tri.file.setGeneNames <- function(filename, geneNames, transcriptNames=NULL){
   trInfo <- tri.load(filename);
   trInfo <- tri.setGeneNames(trInfo, geneNames, transcriptNames);
   tri.save(trInfo, filename);
}

## Save transcript information data from DataFrame into file.
tri.save <- function(trInfo, filename){
   stopifnot( is(trInfo,"DataFrame"));
   if(file.exists(filename)){
      file.remove(filename);
   }
   header <- sprintf("# M %i",dim(trInfo)[1]);
   outF <- file(filename, "w");
   writeLines(header,outF);
   write.table(IRanges::as.data.frame(trInfo), outF, quote=FALSE, row.names=FALSE, col.names=FALSE);
   close(outF);
}
