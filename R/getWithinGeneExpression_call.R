
getWithinGeneExpression <- function(sampleFile, outFile=NULL, trInfoFile=NULL, pretend=FALSE){
   if(is.null(trInfoFile)){
      trInfoFile=.changeExt(sampleFile);
   }
   if(!file.exists(trInfoFile)){
      stop("Please provide valid trInfoFile");
   }
   if(is.null(outFile)){
      fName <- tail(unlist(strsplit(tail( unlist(strsplit(sampleFile,"/")),1),"\\\\")),1)
      if(fName != ""){
         outFile <- tempfile(paste(fName,"-",sep=""),fileext =".WGE");
      }else{
         outFile <- tempfile("dataBS-E-",fileext=".WGE");
      }
   }
   args <- c('getWithinGeneExpression', sampleFile, '--outFile', outFile, '--trInfoFile', trInfoFile);

   if(pretend){
      print(paste(args,collapse=" "))
   }else{
      argc <- length(args);
      result <- .C("_getWithinGeneExpression", as.integer(argc), as.character(args));
   }
   return(outFile);
}
