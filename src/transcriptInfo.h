#ifndef TRANSCRIPTINFO_H
#define TRANSCRIPTINFO_H
#include <string>
#include <vector>

using namespace std;

struct transcriptT {//{{{
   string g,t;
   long l;
   long double effL;
   bool operator< (const transcriptT& d2) const{
      if(g==d2.g)return t<d2.t;
      return g<d2.g;
   }
};//}}}

struct geneT {//{{{
   string name;
   long m;
   vector<long> trs;
};//}}}

class TranscriptInfo{
   private:
      long M,G;
      bool ok,ordered;
      vector<transcriptT> transcripts;
      vector<geneT> genes;
      void setGeneInfo();
   public:
      TranscriptInfo();
      void clearTranscriptInfo();
      TranscriptInfo(string fileName);
      bool readInfo(string fileName);
      bool writeInfo(string fileName, bool force = false);
      bool setInfo(vector<string> gNames,vector<string> tNames, vector<long> lengths);
      bool isOK();
      long getM();
      long getG();
      const vector<long>* getGtrs(long i);
      long L(long i);
      long double effL(long i);
      string trName(long i);
      string geName(long i);
      bool genesOrdered();
      void setEffectiveLength(vector<long double> effL);
      vector<long double> *getShiftedLengths(bool effective = false);
};

#endif
