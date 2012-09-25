
#include "TagAlignments.h"
#include "common.h"

//#define MEM_USAGE

TagAlignments::TagAlignments(){//{{{
   knowNtotal=false;
   knowNreads=false;
   Ntotal=0;
   Nreads=0;
}//}}}
void TagAlignments::init(long Nreads,long Ntotal, long M){//{{{
   currentRead = 0;
   reservedN = 0;
   if(Nreads>0){
      this->Nreads=Nreads;
      knowNreads=true;
      readIndex.reserve(Nreads+2);
   }      
   readIndex.push_back(0);
   
   if(Ntotal>0){
      this->Ntotal=Ntotal;
      knowNtotal=true;
      reservedN = Ntotal+1;
      trIds.reserve(reservedN);
      probs.reserve(reservedN);
   }
   if(M>0){
      this->M=M;
      readsInIsoform.assign(M,-1);
   }else{
      readsInIsoform.clear();
      this->M=0;
   }
}//}}}
void TagAlignments::pushAlignment(long trId, double prob){//{{{
   if(trId>=M){
      M=trId+1;
      readsInIsoform.resize(M,-1);
   }
   if(prob<=0)return;
   if(readsInIsoform[trId] == currentRead){
      // the read has already one alignement to this isoform
     for(long i=readIndex[currentRead];i<(long)trIds.size();i++)
        if(trIds[i] == trId){
           probs[i] += prob;
           break;
        }
   }else{
      if(! knowNtotal){
         // the size of arrays is unknown try to reserve sensible amount of space if we know Nreads
         if(knowNreads && reservedN && ((long)probs.size() == reservedN)){
            // we reached the size of reserved space
            double dens = (double)probs.size() / currentRead; 
            dens *= 1.05; //increase it by 5%
            reservedN =(long)( reservedN + (dens) * (Nreads - currentRead + 1000.0) );
         #ifdef MEM_USAGE
            message("TagAlignments:\n   size: %ld  reserving: %ld  capacity before: %ld\n",probs.size(),reservedN,probs.capacity());
         #endif
            trIds.reserve(reservedN);
            probs.reserve(reservedN);
         #ifdef MEM_USAGE
            message("   capacity after: %ld\n",probs.capacity());
         #endif
         }else if(knowNreads && (! reservedN) && (currentRead == Nreads / 4 )){
            // one quarter in, try to reserve sensible amount of space
            double dens = (double)probs.size() / currentRead; 
            dens *= 1.05; //increase it by 5%
            reservedN =(long)((dens) * (Nreads));
         #ifdef MEM_USAGE
            message("TagAlignments:\n   size: %ld  reserving: %ld  capacity before: %ld\n",probs.size(),reservedN,probs.capacity());
         #endif
            trIds.reserve(reservedN);
            probs.reserve(reservedN);
         #ifdef MEM_USAGE
            message("   capacity after: %ld\n",probs.capacity());
         #endif
         }
      }
      trIds.push_back(trId);
      probs.push_back(prob);
      readsInIsoform[trId] = currentRead;
   }
}//}}}
void TagAlignments::pushRead(){//{{{
   if(readIndex[currentRead] == (int_least32_t) probs.size()){
      // no nonzero alignments for this read
      return;
   }
   currentRead++;
   readIndex.push_back(probs.size());
}//}}}
void TagAlignments::finalizeRead(long &M, long &Nreads, long &Ntotal){//{{{
   M = this->M = readsInIsoform.size();
   Nreads = this->Nreads = readIndex.size() - 1;
   Ntotal = this->Ntotal = probs.size();
#ifdef MEM_USAGE
   message("TagAlignments: readIndex size: %ld  capacity %ld\n",readIndex.size(),readIndex.capacity());
   message("TagAlignments: probs size: %ld  capacity %ld\n",probs.size(),probs.capacity());
#endif
}//}}}
int_least32_t TagAlignments::getTrId(long i) const {//{{{
   if(i<Ntotal)return trIds[i];
   return -1;
}//}}}
double TagAlignments::getProb(long i) const {//{{{
   if(i<Ntotal)return probs[i];
   return 0;
}//}}}
int_least32_t TagAlignments::getReadsI(long i) const {//{{{
   if(i<=Nreads)return readIndex[i];
   return 0;
}//}}}

/*void TagAlignments::setProb(long i,double p){//{{{
   if(i<Ntotal)probs[i]=p;
}//}}}*/
