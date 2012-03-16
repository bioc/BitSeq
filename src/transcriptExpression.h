#ifndef TRANSCRIPTEXPRESSION_H
#define TRANSCRIPTEXPRESSION_H
#include<vector>
#include<string>

using namespace std;

enum TE_FileType{ SAMPLER_MEANS, MEAN_VARIANCE };

struct trExpInfoT{
   long double exp,var;
   long id;
   bool operator< (const trExpInfoT& d2) const{
      return exp<d2.exp;
   }
};

class TranscriptExpression{
   private:
      long M;
      bool logged;
      vector<trExpInfoT> trs;

   public:
      TranscriptExpression();
      TranscriptExpression(string fileName, TE_FileType fileType = SAMPLER_MEANS);
      bool readExpression(string fileName, TE_FileType fileType = SAMPLER_MEANS);
      void doSort(bool reverse = false);
      long getM(){return M;}
      bool isLogged(){return logged;}
      long double exp(long tr){return trs[tr].exp;}
      long double var(long tr){return trs[tr].var;}
      long id(long tr){return trs[tr].id;}
};

#endif
