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
   ## omit last column if its NA
   if(is.na(data[1,dim(data)[2]])){
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

# paste string so that it is within 75ch width and lines are ending with \ (i.e. it's usable in bash)
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
