
getGeneExpression <- function(sampleFile, outFile=NULL, trInfo=NULL, trInfoFile=NULL, pretend=FALSE){
   if(!is.null(trInfo)){
      if(!tri.hasGeneNames(trInfo))stop("The trInfo object does not contain necessary information.");
      trInfoFile <- tempfile("data",fileext=c(".tr"));
      if(!pretend){
         tri.save(trInfo,trInfoFile);
      }else{
         message(paste(c("# Save trInfo object as ",trInfoFile," file.")));
      }
   }else {
      if(is.null(trInfoFile)){
         trInfoFile<-.changeExt(sampleFile);
      }
   }
   if(!file.exists(trInfoFile)){
      stop("Please provide valid transcript information either through trInfo object or trInfoFile file.");
   }
   if(!(pretend || tri.file.hasGeneNames(trInfoFile)))stop(paste(c("File ",trInfoFile," does not contain necessary information.")));
   if(is.null(outFile)){
      fName <- tail(unlist(strsplit(tail( unlist(strsplit(sampleFile,"/")),1),"\\\\")),1)
      if(fName != ""){
         outFile <- tempfile(paste(fName,"-",sep=""),fileext =".GE");
      }else{
         outFile <- tempfile("dataBS-E-",fileext=".GE");
      }
   }
   args <- c('getGeneExpression', sampleFile, '--outFile', outFile, '--trInfoFile', trInfoFile);

   if(pretend){
      writeLines(.specialPaste(args))
   }else{
      argc <- length(args);
      result <- .C("_getGeneExpression", as.integer(argc), as.character(args));
   }
   return(outFile);
}
