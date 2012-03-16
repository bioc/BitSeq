
//using namespace std;

#include "sampler.h"

class CollapsedSampler : public Sampler{
   private:
   vector<long> Z;

   void sampleZ();

   // USING unifromDistribution from class Sampler.
   
   public:

   CollapsedSampler();
   virtual ~CollapsedSampler();
   
   //virtual void init(long n, long m, long samplesTotal, long samplesOut, long Nmap, long Nunmap, const vector<long> &alignI, const vector<TagAlignment> &alignments, const distributionParameters &betaPar, const distributionParameters &dirPar);
   virtual void update();
   virtual void sample();
   
};
