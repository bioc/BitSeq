#include<fstream>
#include<map>

#include "common.h"
#include"transcriptInfo.h"

bool TranscriptInfo::genesOrdered(){//{{{
   return ordered;
}//}}}
bool TranscriptInfo::writeInfo(string fileName, bool force){//{{{
   ofstream trF;
   if(! force){
      ifstream testF(fileName.c_str());
      if(testF.is_open()){
         testF.close();
         return false;
      }
      testF.close();
   }
   trF.open(fileName.c_str(),ios::out | ios::trunc);
   if(! trF.is_open() ) return false;
   trF<<"# M "<<M<<endl;
   for(long i=0;i<M;i++)
      trF<<transcripts[i].g<<" "<<transcripts[i].t<<" "<<transcripts[i].l<<" "<<transcripts[i].effL<<endl;
   trF.close();
   return true;
}//}}}
bool TranscriptInfo::setInfo(vector<string> gNames,vector<string> tNames, vector<long> lengths){//{{{
   if((gNames.size()!=tNames.size())||(tNames.size()!=lengths.size())) return false;
   transcriptT newT;
   M = (long) gNames.size();
   for(long i=0;i<M;i++){
      newT.g=gNames[i];
      newT.t=tNames[i];
      newT.l=lengths[i];
      newT.effL = lengths[i];
      transcripts.push_back(newT);
   }
   setGeneInfo();
   ok = true;
   return ok;
}//}}}
void TranscriptInfo::setGeneInfo(){//{{{
   map<string,long> names;
   geneT tmpG;
   long gi=0,i;
   ordered = true;
   string lastName = "-noname-";
   for(i=0;i<M;i++){
      if(transcripts[i].g == lastName){
         genes[gi].m++;
         genes[gi].trs.push_back(i);
      }else{
         if(names.count(transcripts[i].g) == 0){
            tmpG.name = transcripts[i].g;
            tmpG.m = 1;
            tmpG.trs = vector<long>(1,i);
            genes.push_back(tmpG);
            gi=genes.size()-1;
            names[transcripts[i].g] = gi;
            lastName=transcripts[i].g;
         }else{
            ordered=false;
            //warning("TranscriptInfo: Transcripts of gene %ld are not grouped.\n",transcripts[i].g);
            gi = names[transcripts[i].g];
            genes[gi].m++;
            genes[gi].trs.push_back(i);
         }
      }
   }
   G = genes.size();
   //for(i=0;i<G;i++)message("%s %ld %ld\n",(genes[i].name).c_str(),genes[i].m,genes[i].trs.size());
}//}}}
TranscriptInfo::TranscriptInfo(){ clearTranscriptInfo(); }
void TranscriptInfo::clearTranscriptInfo(){//{{{
   M=G=0;
   ok=false;
   ordered=true;
   transcripts.clear();
   genes.clear();
}//}}}
TranscriptInfo::TranscriptInfo(string fileName){//{{{
   // TranscriptInfo();
   readInfo(fileName);
}//}}}
bool TranscriptInfo::readInfo(string fileName){//{{{
   clearTranscriptInfo();
   ifstream trFile(fileName.c_str());
   if(!trFile.is_open()){
      error("TranscriptInfo: problem reading transcript file.\n");
      return false;
   }
   transcriptT newT;
   while(trFile.good()){
      while(trFile.good() && (trFile.peek()=='#'))
         trFile.ignore(100000000,'\n');
      trFile>>newT.g>>newT.t>>newT.l;
      // read effective length if present
      while((trFile.peek() == '\t')||(trFile.peek() == ' ')) trFile.get();
      if(trFile.peek() == '\n') newT.effL = newT.l;
      else trFile>>newT.effL;
      if(trFile.good())
         transcripts.push_back(newT);
      trFile.ignore(100000000,'\n');
   }
   trFile.close();
   ok = true;
   M = (long)transcripts.size();
   setGeneInfo();
   return ok;
}//}}}
bool TranscriptInfo::isOK(){//{{{
   return ok;
}//}}}
long TranscriptInfo::getM(){//{{{
   return M;
}//}}}
long TranscriptInfo::getG(){//{{{
   return G;
}//}}}
const vector<long>* TranscriptInfo::getGtrs(long i){//{{{
   if(i>G) return NULL;
   return &genes[i].trs;
}//}}}
long double TranscriptInfo::effL(long i){//{{{
   if(ok && (i<M))return transcripts[i].effL;
   return 0;
}//}}}
long TranscriptInfo::L(long i){//{{{
   if(ok && (i<M))return transcripts[i].l;
   return 0;
}//}}}
string TranscriptInfo::trName(long i){//{{{
   if(ok && (i<M))return transcripts[i].t;
   return "";
}//}}}
string TranscriptInfo::geName(long i){//{{{
   if(ok && (i<M))return transcripts[i].g;
   return "";
}//}}}
void TranscriptInfo::setEffectiveLength(vector<long double> effL){//{{{
   if((long)effL.size() != M){
      warning("TranscriptInfo: Wrong array size for effective length adjustment.\n");
      return;
   }
   // Adjust effective length to similar scale as normal length
   long double sumL = 0,sumN = 0,norm;
   for(long i=0;i<M;i++){
      sumN+=effL[i];
      sumL+=transcripts[i].l;
   }
// don't normalize
//   norm = sumL / sumN;
   norm = 1;
   for(long i=0;i<M;i++){
      transcripts[i].effL = effL[i] * norm;
   }
}//}}}
vector<long double> *TranscriptInfo::getShiftedLengths(bool effective){//{{{
   vector<long double> *Ls = new vector<long double>(M+1);
   for(long i=0;i<M;i++){
      if(effective)(*Ls)[i+1] = transcripts[i].effL;
      else (*Ls)[i+1] = transcripts[i].l;
   }
   return Ls;
}//}}}
