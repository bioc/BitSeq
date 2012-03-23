
getGeneExpression <- function(sampleFile, outFile=NULL, trInfoFile=NULL){
   if(is.null(trInfoFile)){
      trInfoFile=.changeExt(sampleFile);
   }
   if(!file.exists(trInfoFile)){
      stop("Please provide valid trInfoFile");
   }
   if(is.null(outFile)){
      fName <- tail(unlist(strsplit(tail( unlist(strsplit(sampleFile,"/")),1),"\\\\")),1)
      if(fName != ""){
         outFile <- tempfile(paste(fName,"-",sep=""),fileext =".GE");
      }else{
         outFile <- tempfile("dataBS-E-",fileext=".GE");
      }
   }
   args <- c('getGeneExpression', sampleFile, '--outFile', outFile, '--trInfoFile', trInfoFile);

   argc <- length(args);
   result <- .C("_getGeneExpression", as.integer(argc), as.character(args));
   return(outFile);
}
