#ifndef TAGALIGNMENT_H
#define TAGALIGNMENT_H

class TagAlignment{
   private:
      long trId;
//      bool strand; // true = forward; false = reverse
      long double prob,lowProb;
   public:
      //TagAlignment(long t=0,bool s=true,long double p = 0,long double lp = 0){
      TagAlignment(long t=0,long double p = 0,long double lp = 0){
         trId=t;
//         strand=s;
         prob=p;
         lowProb = lp;
      }
/*      void setParams(long t=0,bool s=true,long double p = 0,long double lp = 0){
         trId=t;
         strand=s;
         prob=p;
         lowProb = lp;
      }*/
      long getTrId()const {return trId;}
//      bool getStrand()const {return strand;}
      double getProb()const {return prob;}
      void setProb(double p){prob=p;}
      double getLowProb()const {return lowProb;}
}; 

#endif
