#include<sys/time.h>

//using namespace std;

#include "collapsedSampler.h"
#include "common.h"

#define DEBUG(x) 

CollapsedSampler::CollapsedSampler(){ //{{{
//   message("COLLAPSED\n");
}//}}}
CollapsedSampler::~CollapsedSampler(){ //{{{
//   message("COLLAPSED DIE\n");
//   Sampler::~Sampler();
}//}}}
/*void CollapsedSampler::init(long n, long m, long samplesTotal, long samplesOut, long Nmap, long Nunmap, const vector<long> &alignI, const vector<TagAlignment> &alignments, const distributionParameters &betaPar, const distributionParameters &dirPar,long seed){//{{{

   Sampler::init(n,m,samplesTotal,samplesOut,Nmap,Nunmap,alignI,alignments,betaPar,dirPar,long seed);

   Z.assign(n,0);
}//}}}*/
void CollapsedSampler::sampleZ(){//{{{
   DEBUG(message("sampleZ\n");)
   long i,j,k;
   if(Sof(Z)!=Nmap){
      Z.assign(Nmap,0);
      // init Z&C
      for(i=0;i<Nmap;i++){
         //choose random transcript;
         k = (long) (m * uniformDistribution(rng_mt));
         Z[i]=k;
         C[k]++;
      }
   }
#ifdef DoSTATS
   nZ++;
   struct timeval start, end;
   gettimeofday(&start, NULL);
#endif
   vector<double> phi(m,0); 
   // phi of size M should be enough 
   // because of summing the probabilities for each isoform when reading the data
   double probNorm,r,sum,const1;

   const1=m * dir->alpha + Nmap - 1;
   // randomize order: ???
   for(i=0;i<Nmap-1;i++){
      probNorm=0;
      C[Z[i]]--; // use counts without the current one 

      for(j=0, k=(*alignI)[i]; k < (*alignI)[i+1]; j++, k++){
         //message("%ld %lf ",(*alignments)[k].getTrId(),(*alignments)[k].getProb());
         if((*alignments)[k].getTrId() == 0){
            phi[j] = (*alignments)[k].getProb() *
               (beta->beta + Nunmap + C[0]) * 
               (const1 - C[0]); // this comes from division in "false part"
         }else{
            phi[j] = (*alignments)[k].getProb() * 
               (beta->alpha + Nmap - 1 - C[0]) * 
               (dir->alpha + C[ (*alignments)[k].getTrId() ]); 
               /* 
               / (m * dir->alpha + Nmap - 1 - C[0]) ;
               this term was replaced by (const1 - C[0]) 
               and moved into "true part" as multiplication 
               */
         }
         //message("%lf\n",phi[j]);
         probNorm += phi[j];
      }
      r = uniformDistribution(rng_mt);
      r*=probNorm;
      for(j = 0, sum = 0 ; sum<r; j++){
         sum += phi[j];
//         sum += phi[j] / probNorm; instead of each divide do r*probNorm
      }
      if(j==0){
         // e.g. if probNorm == 0
         // assign to noise
         j = (*alignI)[i+1]-(*alignI)[i];
      }
      Z[i] = (*alignments)[ (*alignI)[i] + j - 1 ].getTrId();
      C[ Z[i] ]++;
   }
#ifdef DoSTATS
   gettimeofday(&end, NULL);
   tZ += (end.tv_sec-start.tv_sec)*1000*1000+(end.tv_usec-start.tv_usec);
#endif
   DEBUG(message("Z sampled\n");)
}//}}}

void CollapsedSampler::update(){//{{{
   Sampler::update();

   sampleTheta();

   updateSums();
   if((doLog)&&(save))appendFile();
}//}}}
void CollapsedSampler::sample(){//{{{
   Sampler::sample();

   sampleZ();
}//}}}
