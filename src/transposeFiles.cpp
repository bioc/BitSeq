#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

#include "fileHeader.h"
#include "transposeFiles.h"
#include "common.h"

#define FOR(x,y,n) for(x=y;x<n;x++)
#define FR(x,n) FOR(x,0,n)
#define Sof(x) (long)x.size()

bool transposeFiles(vector<string> inFileNames, string outFileName, bool verbose, string message){
   long M=0,fileN=1,i,j,bufMax,bufN,m,n,totalN,maxN=0,f;
   bool trans=false,transposed=false;
   vector<long> N;
   bufMax=BUFFER_DEFAULT;

   ofstream outFile(outFileName.c_str());
   if(!outFile.is_open()){//{{{
      error("TransposeFile: Unable to open output file\n");
      return 0;
   }//}}}
   //{{{ Opening input
   fileN = Sof(inFileNames);
   ifstream inFile[fileN];
   totalN=0;
   FileHeader fh;
   for(i=0;i<fileN;i++){
      inFile[i].open(inFileNames[i].c_str());
      fh.setFile(&inFile[i]);
      if(!fh.samplesHeader(n,m,trans)){
         error("TransposeFile: Unable to read header of file: %s\n",(inFileNames[i]).c_str());
         return 0;
      }
      if(Sof(N)==0){
         M=m;
         transposed=trans;
         maxN=n;
      }else if((M!=m)||(transposed!=trans)){
         error("TransposeFile: Different number of transcripts or file %s is in wrong format.\n",(inFileNames[i]).c_str());
         return 0;
      }
      outFile<<"# "<<inFileNames[i]<<" "<<n<<endl;
      N.push_back(n);
      if(n>maxN)maxN=n;
      totalN+=n;
   }
   if(bufMax>m)bufMax=m;
   //}}}

   outFile<<message;
   if(!trans)
      outFile<<"# T (M rows,N cols)";
   else 
      outFile<<"# (N rows,M cols)";
   outFile<<"\n# M "<<m<<"\n# N "<<totalN<<endl;
   outFile.precision(9);
   outFile<<scientific;
   if(verbose)message("Transposing files:\n Samples: %ld Transcripts: %ld Buffer size: %ld\n",totalN,m,bufMax);
   if(!trans){ // {{{
      vector< vector<long> > seeks(fileN,vector<long>(maxN,-1));
      vector<vector<string> > valueBuf(bufMax,vector<string>(totalN));
      long lastBuf = 0, done=0;
      bufN=bufMax;
      if(verbose){
         message("Read %ld->%ld",done,bufN+done);
         fflush(stdout);
      }
      FR(f,fileN){
         FR(i,N[f]){
            FR(j,bufN) inFile[f]>>valueBuf[j][lastBuf];
            lastBuf++;
            seeks[f][i]=inFile[f].tellg();
            inFile[f].ignore(10000000,'\n');
         }
      }
      if(verbose)message(" write\n");
      FR(j,bufN){
         FR(i,lastBuf)
            outFile<<valueBuf[j][i]<<" ";
         outFile<<endl;
      }
      lastBuf=0;
      done=bufN;
      while(done<m){
         bufN=bufMax;
         if(m-done<bufMax)bufN=m-done;
         if(verbose){
            message("Read %ld->%ld",done,bufN+done);
            fflush(stdout);
         }
         FR(f,fileN){
            FR(i,N[f]){
               inFile[f].seekg(seeks[f][i]);
               FR(j,bufN) inFile[f]>>valueBuf[j][lastBuf];
               lastBuf++;
               seeks[f][i]=inFile[f].tellg();
            }
         }
         if(verbose)message(" write\n");
         FR(j,bufN){
            FR(i,lastBuf)
               outFile<<valueBuf[j][i]<<" ";
            outFile<<endl;
         }
         lastBuf=0;
         done+=bufN;
      }
      FR(f,fileN)inFile[f].close();
   } // }}}
   else{ // if(trans) {{{
      vector<long> seeks(m,-1);
      vector<vector<string> > valueBuf(m,vector<string>(bufMax));
      long done;
      FR(f,fileN){
         seeks.assign(M,-1);
         done = 0;
         while(done<N[f]){
            bufN=bufMax;
            if(bufN>N[f]-done)bufN=N[f]-done;
            if(verbose){
               message("Read file %ld: %ld->%ld",f,done,bufN+done);
               fflush(stdout);
            }
            FR(j,M){
               if(seeks[j]!=-1)inFile[f].seekg(seeks[j]);
               FR(i,bufN){
                  inFile[f]>>valueBuf[j][i];
               }
               seeks[j]=inFile[f].tellg();
               if((j+1<M)&&(seeks[j+1]==-1))inFile[f].ignore(100000000,'\n');
            }
            if(verbose)message(" write\n");
            FR(i,bufN){
               FR(j,M)
                  outFile<<valueBuf[j][i]<<" ";
               outFile<<endl;
            }
            done+=bufN;
         }
         inFile[f].close();      
      }
   } //}}}
   outFile.close();
   return 1;
}
