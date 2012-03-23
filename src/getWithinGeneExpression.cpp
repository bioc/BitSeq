/*
 *
 * Produce relative expression within gene
 *
 */
#include<cmath>
#include<cstdio>

using namespace std;

#include "posteriorSamples.h"
#include "transcriptInfo.h"
#include "common.h"
#include "argumentParser.h"

#define Sof(x) (long)x.size()
#define FOR(x,y,n) for(x=y;x<n;x++)
#define FR(x,n) FOR(x,0,n)

extern "C" int getWithinGeneExpression(int *argc,char* argv[]){
   buildTime(argv[0]);
   string programDescription=
"Computes relative expression of transcripts within genes.\n\
   [sampleFiles] should contain transposed MCMC.\n\
   program can produce means and variance and write them into [sumFile]\n\
   or individual MCMC samples which are written into [outFile].";   
   // Set options {{{
   ArgumentParser args(programDescription,"[sampleFiles]",1);
   args.addOptionS("t","trInfoFile","trInfoFileName",1,"Name of the transcript file.");
   args.addOptionB("a","adjustByLength","adjust",0,"Adjust expression by transcripts length.");
   args.addOptionS("o","outFile","outFileName",0,"Name of the output file.");
   args.addOptionS("s","sumFile","sumFileName",0,"Name of summarization file where true mean, true variance and relative mean and realtive variance are saved.");
   args.addOptionB("l","log","log",0,"Use logged values.");
//   args.addOptionS("t","type","type",0,"Type of variance, possible values: [sample,sqDif] for sample variance or sqared difference.","sample");
   if(!args.parse(*argc,argv))return 0;
   // }}}
   bool doLog=args.flag("log"),doOut=args.isSet("outFileName"),doSum=args.isSet("sumFileName");
   if(! (doOut || doSum)){
      error("Main: Have to secify at least one of -o -s.\n");
      return 1;
   }

   if(doLog){
      if(args.verbose)message("Using logged values.\n");
   }else{
      if(args.verbose)message("NOT using logged values.\n");
   }
   
   long i,j,N=0,M=0,G,g,gM,m;
   TranscriptInfo trFile;
   if(! trFile.readInfo(args.getS("trInfoFileName"))){
      error("Main: Failed loading transcript file.\n");
      return 1;
   }
   G = trFile.getG();

   PosteriorSamples  samples;
   if((! samples.initSet(M,N,args.args()[0] ))||(M<=0)||(N<=0)){
      error("Main: Failed loading MCMC samples.\n");
      return 1;
   }
   if(M!=trFile.getM()){
      error("Main: Number of transcripts does not match: %ld %ld\n",trFile.getM(),M);
      return 1;
   }
   if(args.verbose)message("Genes: %ld\nTranscripts: %ld\n",G,M);

   ofstream outFile,sumFile;
   if(doOut){
      outFile.open(args.getS("outFileName").c_str());
      if(! outFile.is_open()){
         error("Main: %s file write failed!\n",(args.getS("outFileName")).c_str());
         return 1;
      }
      outFile<<"# from: "<<args.args()[0]<<"\n# samples of within gene expression\n";
      if(! trFile.genesOrdered()){
         message("WARNING: transcripts in output file will be reordered and grouped by genes.\n");
         outFile<<"# WARNING: transcripts in output file are be reordered and grouped by genes.\n";
      }
      if(doLog)outFile<<"# L \n";
      outFile<<"# T (M rows,N cols)\n";
      outFile<<"# M "<<M<<"\n# N "<<N<<endl;
   }
   if(doSum){
      sumFile.open(args.getS("sumFileName").c_str());
      if(! sumFile.is_open()){
         error("Main: %s file write failed!\n",(args.getS("sumFileName")).c_str());
         return 1;
      }
      sumFile<<"# from: "<<args.args()[0]<<"\n# <mean> <variance> <mean of within gene expression>  <variance of within gene expression>\n# M "<<M<<endl;
   }
   vector<long double> mean(M),mean2(M),sqSum(M),sqSum2(M);
   vector< vector<double> > trs;
   vector<long double> normals(N,0);
   long double x,sum,var,var2,l;
   if(args.flag("adjust")&&(doSum)){
      vector<double> tr(M);
      if(args.verbose)message("Computing normalization constants, because of length adjustment.\n");
      FR(j,M){
         if(args.verbose)progressLog(j,M);
         samples.getTranscript(j,tr);
         FR(i,N)
            normals[i] += tr[i]/trFile.L(j);
      }
   }
   if(args.verbose)message("Computing within gene relative expression.\n");
   FR(g,G){
      if(args.verbose)progressLog(g,G);
      gM = trFile.getGtrs(g)->size();
      if(Sof(trs)<gM)trs.resize(gM);
      //message("%ld\n",gM);
      FR(j,gM){
         m = (*trFile.getGtrs(g))[j];
         samples.getTranscript( m , trs[j]);
         mean[m] = mean2[m] = sqSum[m] = sqSum2[m] = 0;
      }
      FR(i,N){
         sum = 0;
         FR(j,gM){
            if(args.flag("adjust")){
               m = (*trFile.getGtrs(g))[j];
               l = trFile.L(m);
               trs[j][i]/=l;
               if(normals[i]!=0)trs[j][i]/=normals[i];
            }
            x=trs[j][i];
            sum += x;
            if(doSum){
               m = (*trFile.getGtrs(g))[j];
               if(doLog)x=log(trs[j][i]);
               mean[m] += x;
               sqSum[m] += x*x;
            }
         }
         FR(j,gM){
            trs[j][i] /= sum;
            if(doLog)trs[j][i] = log(trs[j][i]);
            if(doSum){
               m = (*trFile.getGtrs(g))[j];
               mean2[m] += trs[j][i];
               sqSum2[m] += trs[j][i]*trs[j][i];
            }
         }
      }
      if(doOut){
         FR(j,gM){
            FR(i,N)
               outFile<<trs[j][i]<<" ";
            outFile<<endl;
         }
      }
   }
   if(doOut)outFile.close();
   if(doSum){
      FR(j,M){
         mean[j] /= N;
         var = sqSum[j]/N - mean[j]*mean[j];
         mean2[j] /= N;
         var2 = sqSum2[j]/N - mean2[j]*mean2[j];
         sumFile<<mean[j]<<" "<<var<<" "<<mean2[j]<<" "<<var2<<endl;
      }
      sumFile.close();
   }
   if(args.verbose)message("DONE\n");
   return 0;
}


#ifndef BIOC_BUILD
int main(int argc,char* argv[]){
   return getWithinGeneExpression(&argc,argv);
}
#endif
