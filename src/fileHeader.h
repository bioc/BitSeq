#include<cstdio>
#include<fstream>
#include<sstream>
#include<string>

using namespace std;

#include "common.h"

class FileHeader{
   private:
      ifstream *file;
   public:
      FileHeader(ifstream *f=NULL){//{{{
         file=f;
      }//}}}
      void setFile(ifstream *f){//{{{
         file=f;
      }//}}}
      bool samplesHeader(long &n,long &m,bool &transposed,bool &logged){//{{{
         transposed=false;
         logged = false;
         if((file==NULL)||(!file->is_open())){
            error("No file for header read.\n");
            m=0;
            n=0;
            return false;
         }
         string line,str;
         istringstream lineS;
         while((!file->eof())&&(file->peek() == '#')){
            getline(*file, line);
            while((!file->eof())&&((file->peek() == ' ')||(file->peek() == '\n')))file->get();
            lineS.clear();
            lineS.str(line);
            while(lineS.good()){
               lineS>>str;
               if(str == "M"){
                  lineS>>m;
//                  printf("header %ld\n",m);
                  continue;
               }
               if(str == "N"){
                  lineS>>n;
//                  printf("header %ld\n",n);
                  continue;
               }
               if(str == "T"){
//                  printf("header transposed\n");
                  transposed = true;
                  continue;
               }  
               if(str == "L"){
//                  printf("header transposed\n");
                  logged = true;
                  continue;
               }  
            }
         }
         return true;
      }//}}}
      bool samplesHeader(long &n,long &m,bool &transposed){//{{{
         transposed=false;
         if((file==NULL)||(!file->is_open())){
            error("No file for header read.\n");
            m=0;
            n=0;
            return false;
         }
         string line,str;
         istringstream lineS;
         while((!file->eof())&&(file->peek() == '#')){
            getline(*file, line);
            while((!file->eof())&&((file->peek() == ' ')||(file->peek() == '\n')))file->get();
            lineS.clear();
            lineS.str(line);
            while(lineS.good()){
               lineS>>str;
               if(str == "M"){
                  lineS>>m;
//                  printf("header %ld\n",m);
                  continue;
               }
               if(str == "N"){
                  lineS>>n;
//                  printf("header %ld\n",n);
                  continue;
               }
               if(str == "T"){
//                  printf("header transposed\n");
                  transposed = true;
                  continue;
               }  
            }
         }
         return true;
      }//}}}
      bool transcriptsHeader(long &m){//{{{
         if((file==NULL)||(!file->is_open())){
            error("No file for header read.\n");
            m=0;
            return false;
         }
         string line,str;
         istringstream lineS;
         while((!file->eof())&&(file->peek() == '#')){
            getline(*file, line);
            while((!file->eof())&&((file->peek() == ' ')||(file->peek() == '\n')))file->get();
            lineS.clear();
            lineS.str(line);
            while(lineS.good()){
               lineS>>str;
               if(str == "M"){
                  lineS>>m;
//                  printf("header %ld\n",m);
                  continue;
               }
            }
         }
         return true;
      }//}}}
      bool transcriptsHeader(long &m,long &colN){//{{{
         if((file==NULL)||(!file->is_open())){
            error("No file for header read.\n");
            m=0;
            return false;
         }
         string line,str;
         istringstream lineS;
         while((!file->eof())&&(file->peek() == '#')){
            getline(*file, line);
            while((!file->eof())&&((file->peek() == ' ')||(file->peek() == '\n')))file->get();
            lineS.clear();
            lineS.str(line);
            while(lineS.good()){
               lineS>>str;
               if(str == "M"){
                  lineS>>m;
//                  printf("header %ld\n",m);
                  continue;
               }
               if(str == "colN"){
                  lineS>>colN;
//                  printf("header %ld\n",m);
                  continue;
               }
            }
         }
         return true;
      }//}}}
      bool probHeader(long &Nmap,long &Ntotal,bool &newformat){//{{{
         if((file==NULL)||(!file->is_open())){
            error("No file for header read.\n");
            Nmap=0;
            return false;
         }
         newformat = false;
         string line,str;
         istringstream lineS;
         while((!file->eof())&&(file->peek() == '#')){
            getline(*file, line);
            while((!file->eof())&&((file->peek() == ' ')||(file->peek() == '\n')))file->get();
            lineS.clear();
            lineS.str(line);
            while(lineS.good()){
               lineS>>str;
               if(str == "Ntotal"){
                  lineS>>Ntotal;
//                  printf("header %ld\n",m);
                  continue;
               }
               if(str == "Nmap"){
                  lineS>>Nmap;
//                  printf("header %ld\n",m);
                  continue;
               }
               if(str == "NEWFORMAT"){
                  newformat = true;
               }
            }
         }
         return true;
      }//}}}
      bool varianceHeader(long &m,bool &logged){//{{{
         if((file==NULL)||(!file->is_open())){
            error("No file for header read.\n");
            m=0;
            return false;
         }
         string line,str;
         istringstream lineS;
         while((!file->eof())&&(file->peek() == '#')){
            getline(*file, line);
            while((!file->eof())&&((file->peek() == ' ')||(file->peek() == '\n')))file->get();
            lineS.clear();
            lineS.str(line);
            while(lineS.good()){
               lineS>>str;
               if(str == "M"){
                  lineS>>m;
//                  printf("header %ld\n",m);
                  continue;
               }
               if(str == "L"){
//                  printf("header transposed\n");
                  logged = true;
                  continue;
               }  
            }
         }
         return true;
      }//}}}
      bool paramsHeader(long &parN){//{{{
         if((file==NULL)||(!file->is_open())){
            error("No file for header read.\n");
            parN=0;
            return false;
         }
         string line,str;
         istringstream lineS;
         parN = 0;
         while((!file->eof())&&(file->peek() == '#')){
            getline(*file, line);
            while((!file->eof())&&((file->peek() == ' ')||(file->peek() == '\n')))file->get();
            lineS.clear();
            lineS.str(line);
            while(lineS.good()){
               lineS>>str;
               if(str == "PN"){
                  lineS>>parN;
//                  printf("header %ld\n",parN);
                  continue;
               }
            }
         }
         return true;
      }//}}}
      void close(){//{{{
         file->close();
         file=NULL;
      }//}}}
};
