#include <cstdlib>
#include <string>

#include "common.h"

using namespace std;

void buildTime(char *argv0){
   return ; // dont want to print compile information
   string compileDate = __DATE__;
   string compileTime = __TIME__;
   message("### %s build: %s %s\n",argv0,compileDate.c_str(),compileTime.c_str());
}

bool progressLog(long cur,long outOf, long steps) {
   // output progress status every 10%
   if((outOf>steps)&&(cur%((long)(outOf/steps))==0)&&(cur!=0)){
      message("# %ld done.",cur);
      return true;
   }
   return false;
}
