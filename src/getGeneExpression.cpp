/*
 *
 * Produce overall gene expression
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

extern "C" int getGeneExpression(int *argc,char* argv[]){
   buildTime(argv[0]);
   string programDescription=
"Computes expression of whole genes.\n\
   [sampleFiles] should contain transposed MCMC samples which are transformed into gene expression samples.";   
   // Set options {{{
   ArgumentParser args(programDescription,"[sampleFiles]",1);
   args.addOptionS("t","trInfoFile","trInfoFileName",1,"Name of the transcript file.");
   args.addOptionB("a","adjustByLength","adjust",0,"Adjust expression by transcripts length.");
   args.addOptionB("","rpkm","rpkm",0,"Output RPKM.");
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   args.addOptionB("l","log","log",0,"Output logged values.");
   if(!args.parse(*argc,argv))return 0;
   // }}}
   bool doLog=args.flag("log"),doAdjust=args.flag("adjust")||args.flag("rpkm"),doRPKM=args.flag("rpkm");

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
   if((! samples.initSet(M,N,args.args()[0]))||(M<=0)||(N<=0)){
      error("Main: Failed loading MCMC samples.\n");
      return 1;
   }
   if(M!=trFile.getM()){
      error("Main: Number of transcripts does not match: %ld %ld\n",trFile.getM(),M);
      return 1;
   }
   if(args.verbose)message("Genes: %ld\nTranscripts: %ld\n",G,M);

   ofstream outFile;
   outFile.open(args.getS("outFileName").c_str());
   if(! outFile.is_open()){
      error("Main: %s file write failed!\n",(args.getS("outFileName")).c_str());
      return 1;
   }
   outFile<<"# from: "<<args.args()[0]<<"\n# samples of gene expression\n";
   if(args.verbose)message("Genes will be ordered as they first appear in %s.\n",(args.getS("trInfoFileName")).c_str());
   outFile<<"# Genes will be ordered as they first appear in "<<args.getS("trInfoFileName")<<"\n";
   if(doRPKM)outFile<<"# data in RPKM\n";
   if(doLog)outFile<<"# L \n";
   outFile<<"# T (M rows,N cols)\n";
   outFile<<"# G = M "<<G<<"\n# N "<<N<<endl;
   vector< vector<double> > trs;
   vector<long double> normals(N,0);
   long double sum;
   if(doAdjust){
      vector<double> tr(M);
      if(args.verbose)message("Computing normalization constants, because of length adjustment.\n");
      FR(j,M){
         if(args.verbose)progressLog(j,M);
         samples.getTranscript(j,tr);
         FR(i,N)
            normals[i] += tr[i]/trFile.L(j);
      }
   }
   if(args.verbose)message("Computing gene expression.\n");
   FR(g,G){
      if(args.verbose)progressLog(g,G);
      gM = trFile.getGtrs(g)->size();
      if(Sof(trs)<gM)trs.resize(gM);
      //message("%ld\n",gM);
      FR(j,gM){
         m = (*trFile.getGtrs(g))[j];
         samples.getTranscript( m , trs[j]);
      }
      FR(i,N){
         sum = 0;
         FR(j,gM){
            if(doAdjust&&(normals[i]>0)){
               m = (*trFile.getGtrs(g))[j];
               sum+=(trs[j][i] / trFile.L(m)) / normals[i];
            }else
               sum+=trs[j][i];
         }
         if(doRPKM)sum=sum*10e9;
         if(doLog)sum=log(sum);
         outFile<<sum<<" ";
      }
      outFile<<endl;
   }
   outFile.close();
   if(args.verbose)message("DONE\n");
   return 0;
}

#ifndef BIOC_BUILD
int main(int argc,char* argv[]){
   return getGeneExpression(&argc,argv);
}
#endif
