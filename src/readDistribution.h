#ifndef READDISTRIBUTION_H
#define READDISTRIBUTION_H

#include<vector>
#include<map>

using namespace std;

#include "transcriptInfo.h"
#include "transcriptSequence.h"
#include "transcriptExpression.h"

//#include "samtools/bam.h"
#include "samtoolsHeader.h"

#ifdef BIOC_BUILD
#include <Rinternals.h>
#define bam_init1() ((bam1_t*)S_alloc(1, sizeof(bam1_t)))
// empty destroy, R frees memory itself
#define bam_destroy1(b) 

#else
#define bam_init1() ((bam1_t*)calloc(1, sizeof(bam1_t)))
#define bam_destroy1(b) do { \
   if (b) { free((b)->data); free(b); }	\
} while (0)

#endif


// Defaults: {{{
#define LOW_PROB_MISSES 6
#define MAX_NODE_PAR 2
const long trSizes [] = { 1334,2104,2977,4389};
const long trSizesN = 4;
const long trNumberOfBins = 20;
const long vlmmNodeDependence [] = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0};
const long vlmmNodesN = 21;
const long vlmmStartOffset = 8;
const long pows4 [] = {1,4,16,64,256,1024,4096};
//}}}

struct fragmentT{//{{{
   bam1_t *first,*second;
   bool paired;
   fragmentT(){
      first = bam_init1();
      second = bam_init1();
      paired = true;
   }
   ~fragmentT(){
      bam_destroy1(first);
      bam_destroy1(second);
   }
};

typedef fragmentT* fragmentP;
//}}}

class VlmmNode{//{{{
   private:
      long parentsN;
      vector<long double> probs;

   public:
      VlmmNode(){parentsN = 0;}
      VlmmNode(long p);
      void setParentsN(long p);
      void update(long double Iexp, char b, char bp, char bpp);
      void normalize();
      long double getP(char b, char bp, char bpp);
      long double getPsum(char b);
};//}}}

enum biasT { readM_5, readM_3, uniformM_5, uniformM_3, weight_5, weight_3};
enum readT { mate_5, mate_3, FullPair };

class ReadDistribution{
   private:
      long M,fragSeen,singleReadLength,minFragLen;
      long double lMu,lSigma,logLengthSum,logLengthSqSum;
      long lowProbMismatches;
      bool verbose,uniform,lengthSet,gotExpression,normalized;
      bool warnPos, validLength;
      TranscriptInfo* trInf;
      TranscriptSequence* trSeq;
      TranscriptExpression* trExp;
      // for each transcript, remember seen fragments:
      vector<map<long,long double> > trFragSeen5,trFragSeen3;
      // cache for already computed weight norms for single reads 4',3', Pair x Transcript x Length
      vector<vector<map<long, long double> > > weightNorms;
      // position probability arrays (RE-FACTOR to array of 4 vectors)
      vector<vector<vector<long double> > > posProb;
      vector<vector<VlmmNode> > seqProb;
   
      long double getLengthP(long double len);
      long double getLengthNorm(long double trLen);
      void updatePosBias(long pos, biasT bias, long tid, long double Iexp);
      void updateSeqBias(long pos, biasT bias, long tid, long double Iexp);
      long double getPosBias(long pos, readT read, long tid);
      long double getSeqBias(long pos, readT read, long tid);
      long double getWeightNorm(long len, readT, long tid);
      pair<long double, long double> getSequenceProb(bam1_t *samA);
   public:
      ReadDistribution(long m);
      bool init(TranscriptInfo* trI, TranscriptSequence* trS, TranscriptExpression* trE, bool verb = true);
      bool initUniform(TranscriptInfo* trI, TranscriptSequence* trS, bool verb = true);
      void setLowProbMismatches(long m);
      void setLength(long double mu, long double sigma);
      void observed(fragmentP frag);
      void normalize();
      void logProfiles(string logFileName = "");
      void getP(fragmentP frag,long double &prob,long double &probNoise);
      vector<long double> getEffectiveLengths();
}; 

#endif
