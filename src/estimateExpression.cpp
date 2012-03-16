#include<omp.h>
#include<cmath>
#include<algorithm>

//#include "gibbsParameters.h"
#include "collapsedSampler.h"
#include "gibbsSampler.h"
#include "sampler.h"
#include "fileHeader.h"
#include "myTimer.h"
#include "argumentParser.h"
#include "transcriptInfo.h"
#include "transposeFiles.h"
#include "common.h"

//#define Sof(x) (long)x.size()

//using namespace std;

#define DEBUG(x)
#define FF first
#define SS second

vector<TagAlignment> alignments;
vector<long> alignI; // index of read's first hit in hits
TranscriptInfo trInfo;

long  M;//, mAll; // M : number of transcripts (include transcript 0 ~ Noise)
long N, Nunmap, Nmap; // N: number of read, unmappable read, mappable reads
//double ratioNSoverT;

string outTypeS;
outputType outTypeI;
vector<string> samplesFileNames;
string failedMessage;

void clearDataEE(){
   alignments.clear();
   alignI.clear();
   samplesFileNames.clear();
}

void readData(ArgumentParser &args) {//{{{
   long i,j,k,num,tid;
   double prb;
   long Ntotal=0;
   string readId,strand,blank;
   ifstream inFile;
   MyTimer timer;

   // Read alignment probabilities {{{
   inFile.open(args.args()[0].c_str());
   FileHeader fh(&inFile);
   bool newformat=true;
   if((!fh.probHeader(Nmap,Ntotal,newformat)) || (Nmap ==0)){//{{{
     error("Prob file header read failed.\n");
   }//}}}
   message("Mappings: %ld\n",Nmap);
   message("Ntotal: %ld\n",Ntotal);
   if(Ntotal>Nmap)Nunmap=Ntotal-Nmap;
   else Nunmap=1; //no valid count file assume only one not aligned properly
   N=Nmap+Nunmap;
   alignI.push_back(0);
   vector<long> readInIsoform;
   if(M>0){
      readInIsoform.resize(M,-1);
   }
   long mod=10;
   long isoformsHit=0;
   //long bad = 0;
   timer.start();
   for(i = 0; i < Nmap; i++) {
      inFile>>readId>>num;
      if(!newformat)inFile>>blank;
      if(!inFile.good())break;
     //    message("%s %ld\n",(readId).c_str(),num);
      for(j = 0; j < num; j++) {
         if(newformat) inFile>>tid>>prb;
         else inFile>>tid>>strand>>prb;
         if(inFile.fail()){
            inFile.clear();
            // ignore rest of line
            //inFile.ignore(10000000,'\n');
            // ignore other read's alignments
            j=num;
            // this read goes to noise
            tid=0;
            prb=100;
            //bad++;
         }
         if((!newformat) && (tid!=0)){
            // these lengths are not shifted
            prb /= trInfo.L(tid-1);
         }
         
         // for cases when M is initially unknown
         if(tid>M){
            M=tid+1;
            readInIsoform.resize(M,-1);
         }
         if( readInIsoform[tid] == i){
            // this means that this read already had hit in isoform TID
            // we find this hit and add probability, don't care about the direction
            for(k = Sof(alignments)-1; alignments[k].getTrId() != tid; k--);
            //assert(k >= alignI[i]);
            alignments[k].setProb(alignments[k].getProb() + prb);
         }else{
            if(readInIsoform[tid] == -1)isoformsHit++;
            // note that this read has hit in isoform SID
            readInIsoform[tid] = i;
            alignments.push_back(TagAlignment( tid, prb));
         }
      }
      // ignore rest of line
      inFile.ignore(10000000,'\n');
      alignI.push_back(Sof(alignments)); 
      //assert(alignI[i + 1] - alignI[i] >= 2);

      R_CheckUserInterrupt();
      if((i % mod == 0)&&(i>0)){
         message("  %ld ",i);
         timer.split();
         mod*=10;
      }
   }
   inFile.close();
   //}}}
   if(i<Nmap)message("Read only %ld mappings.\n",i);
   message("Finished Reading!\nTotal hits = %ld\n",alignI[Nmap]);
   message("Isoforms hit: %ld (out of %ld)\n",isoformsHit,M);
   /* {{{ remapping isoforms to ignore those without any hits
   M = mAll;
   M = isoformsHit;
   isoformMap.assign(M);
   for(i=0,j=0;i<mAll;i++)
      if(readInIsoform[i]!=-1){
         readInIsoform[i]=j;
         isoformMap[j]=i;
         j++;
      }
   for(i=0;i<Sof(alignments);i++){
      alignments[i].setTrId( readInIsoform[ alignments[i].getTrId() ] );
   }
   }}}*/
}//}}}

ofstream rLog;

void MCMC(gibbsParameters &gPar,ArgumentParser &args){//{{{
   // Declarations: {{{
   DEBUG(message("Declaratinos:\n"));
   long i,j,samplesHave=0,totalSamples=0,samplesN,chainsN,samplesSave;
   pairD rMean,tmpA,tmpV;
   double rH1,rH2;
   ofstream meansFile,samplesFile[gPar.chainsN()];
   MyTimer timer;
   bool quitNext = false;
   vector<pairD> betwVar(M),withVar(M),s2j(M),totAverage(M),av,var;
   vector<pair<pairD,long> > rHat2(M);
   // }}}
   // Init: {{{
   DEBUG(message("Initialization:\n"));
   samplesN=gPar.samplesN();
   chainsN=gPar.chainsN();
   samplesSave=gPar.samplesSave();

   vector<Sampler*> samplers(chainsN);
   if( ! args.flag("gibbs")){
      for(j=0;j<chainsN;j++)
         samplers[j] = new CollapsedSampler;
   }else{
      for(j=0;j<chainsN;j++)
         samplers[j] = new GibbsSampler;
   }

   // don't want to create a new one every time parallel code starts

   timer.start();;
   // parallel block: 
   // make sure that all functions used are CONST and variables are beaing READ or private
   // private: samplesHave
   #pragma omp parallel for private(samplesHave)
   for(i=0;i<chainsN;i++){
      // Init samplers
      DEBUG(message("Sampler %ld init and burn in.\n",i);)
      samplers[i]->noSave();
      DEBUG(message("init\n");)
      samplers[i]->init(M, samplesN, samplesSave, Nmap, Nunmap,  alignI, alignments, gPar.beta(), gPar.dir(), i*time(NULL));

      DEBUG(message(" burn in\n");) 
      for(samplesHave=0;samplesHave<gPar.burnIn();samplesHave++){
        samplers[i]->sample();
	R_CheckUserInterrupt();
      }
   } 
   message("Burn in: %ld DONE. ",gPar.burnIn());
   timer.split(0,'m');
   //}}}
   // Main sampling loop:
   while(1){
      timer.start();
      // Sample: {{{
      // parallel block:
      // make sure that all functions used are CONST and variables are being READ or private
      // private: samplesHave, samplesSkipped
      #pragma omp parallel for private(samplesHave)
      for(i=0;i<chainsN;i++){
         for(samplesHave = 0;samplesHave<samplesN;samplesHave++){
            samplers[i]->sample();
            samplers[i]->update();
	    R_CheckUserInterrupt();
         }
      }
      totalSamples+=samplesN*chainsN;
      message("\nSampling DONE. ");
      //timeWas=timer.split(0,'m'); only for rLog
      timer.split(0,'m');
      //}}}
      // Check for change of parameters: {{{
      gPar.readParameters();
      // }}}
      // Compute convergence statistics {{{
      totAverage.assign(M,pairD(0,0));
      betwVar.assign(M,pairD(0,0));
      withVar.assign(M,pairD(0,0));
      samplesHave = samplesN;
      for(i=0;i<M;i++){
         for(j=0;j<chainsN;j++){
            tmpA = samplers[j]->getAverage(i);
            tmpV = samplers[j]->getWithinVariance(i);
            totAverage[i].FF += tmpA.FF;
            totAverage[i].SS += tmpA.SS;
            withVar[i].FF += tmpV.FF;
            withVar[i].SS += tmpV.SS;
         }
         totAverage[i].FF /= chainsN;
         totAverage[i].SS /= chainsN;
         withVar[i].FF /= chainsN;
         withVar[i].SS /= chainsN;
         for(j=0;j<chainsN;j++){
            tmpA = samplers[j]->getAverage(i);
            betwVar[i].FF += (totAverage[i].FF - tmpA.FF)*(totAverage[i].FF - tmpA.FF);
            betwVar[i].SS += (totAverage[i].SS - tmpA.SS)*(totAverage[i].SS - tmpA.SS);
         }
         betwVar[i].FF /= (chainsN-1.0);
         betwVar[i].SS /= (chainsN-1.0);
      }
      for(i=0;i<M;i++){
         // betwVar[i] *= samplesHave / (chainsN - 1.0);
         rHat2[i].SS=i;
         if(withVar[i].FF == 0 ){
            rHat2[i].FF.FF = 0;
            rHat2[i].FF.SS = 0;
         } else { 
            rHat2[i].FF.FF = (samplesHave - 1.0) / samplesHave + betwVar[i].FF / withVar[i].FF ;
            rHat2[i].FF.SS = (samplesHave - 1.0) / samplesHave + betwVar[i].SS / withVar[i].SS ;
               //betwVar[i] / ( samplesHave * withVar[i] );
         }
      }
      sort(rHat2.rbegin(),rHat2.rend());
      message("rHat (for %ld samples) \n",samplesN);
      rMean.FF=0;
      rMean.SS=0;
      message("   rHat (rHat from subset | tid | mean theta)\n");
      for(i=0;(i<10) && (i<M);i++){
         rH1 = sqrt(rHat2[i].FF.FF);
         rH2 = sqrt(rHat2[i].FF.SS);
         rMean.FF+=rH1;
         rMean.SS+=rH2;
//         message("   %lf (%lf | %ld | %lf|%lf|%lf)",rHat2[i].FF.FF,rHat2[i].FF.SS,rHat2[i].SS,totAverage[rHat2[i].SS].FF,withVar[rHat2[i].SS].FF,betwVar[rHat2[i].SS].FF/samplesHave);
         message("   %lf (%lf | %ld | %lf)",rH1,rH2,rHat2[i].SS-1,totAverage[rHat2[i].SS].FF);
         message("\n");
//                  message("   %lf",sqrt(rHat2[i].FF));
      }                  
      rMean.FF /= 10.0;
      rMean.SS /= 10.0;
      message("  Mean rHat of worst 10 transcripts.(target: %lf)\n",gPar.targetScaleReduction());
      message("   %lf\n",rMean.FF);
      message("  Mean C0: (");
      for(j=0;j<chainsN;j++)message("%ld ",samplers[j]->getAverageC0());
      message("). Nunmap: %ld\n",Nunmap);
      if(args.flag("gibbs"))message("  Mean thetaAct (noise parameter)\n   %lf\n",totAverage[0].FF);
      message("\n");
//      rLog<<rMean.FF<<" "<<totalSamples<<" "<<timeWas<<" "<<rMean.SS<<endl;
      //}}}
      // Increase sample size and start over: {{{
      if(quitNext){
         if(sqrt(rHat2[0].FF.FF) > gPar.targetScaleReduction()){
            message("WARNING: Following transcripts failed to converge entirely\n   (however the estimates might still be usable):\n");
            stringstream sstr;
            sstr.str();
            sstr<<"# unconverged_transcripts: ";
            for(i=0;(i<M) && (sqrt(rHat2[i].FF.FF) > gPar.targetScaleReduction());i++){
               sstr<<rHat2[i].SS<<" ("<<sqrt(rHat2[i].FF.FF)<<") ";
               message("   %s( %ld , %lf )\n",(trInfo.trName(rHat2[i].SS-1)).c_str(),rHat2[i].SS-1,sqrt(rHat2[i].FF.FF));
            }
            sstr<<"\n";
            failedMessage=sstr.str();
         }
         for(j=0;j<chainsN;j++){
            samplers[j]->noSave();
            samplesFile[j].close();
         }
         break;
      }
      if( (rMean.FF < gPar.targetScaleReduction())||
          (samplesN*2>=gPar.samplesNmax())||
          (totalSamples>5000000) ){
         // Prepare for producing samples if Rhat^2<target scale reduction
         // OR reached samplesNmax
         // OR produced too many samples (>500 000)
         if(! ((rMean.FF < gPar.targetScaleReduction())||
               (totalSamples>5000000)) ){
            samplesN = gPar.samplesNmax();
         }
         quitNext = true;
         message("Producing %ld final samples.\n",samplesN);
         stringstream sstr;
         for(j=0;j<chainsN;j++){
            R_CheckUserInterrupt();
            sstr.str("");
            sstr<<args.getS("outFilePrefix")<<"."<<outTypeS<<"S-"<<j;
            samplesFileNames.push_back(sstr.str());
            samplesFile[j].open(samplesFileNames[j].c_str());
            if(! samplesFile[j].is_open()){
               error("Main: Unable to open output file '%s'.\n",(sstr.str()).c_str());
            }else{
               samplesFile[j]<<"#\n# M "<<M-1<<"\n# N "<<min(samplesSave,samplesN)<<endl;
               samplers[j]->saveSamples(&samplesFile[j],trInfo.getShiftedLengths(true),outTypeI);
            }
         }
      }else{
         samplesN *= 2;
         if(samplesN>gPar.samplesNmax()){
            samplesN = gPar.samplesNmax();
         }
      }
      for(j=0;j<chainsN;j++){
         samplers[j]->resetSampler(samplesN);
      }
      samplesHave=0;
      //}}}
   }
   // Write means: {{{
   meansFile.open((args.getS("outFilePrefix")+".thetaMeans").c_str());
   meansFile<<"# T => Mrows Ncols\n# M "<<M-1<<"\n# N "<<chainsN+1<<endl;
   meansFile<<"# thetaAct variance:"; //prefix 0 isoform with #
   meansFile<<scientific;
   meansFile.precision(9);
   double sum,sum2;
   for(i=0;i<M;i++){
      sum=0;
      sum2=0;
      for(j=0;j<chainsN;j++){
         sum+=samplers[j]->getAverage(i).FF; 
         sum2+=samplers[j]->getAverage(i).SS; 
         if(i==0)meansFile<<" "<<samplers[j]->getWithinVariance(0).FF;
      }
      if(i==0){
         meansFile<<"\n# thetaAct "<<sum/chainsN<<" ";
      }else{
         meansFile<< i<<" "<< sum/chainsN <<" ";
      }
      for(j=0;j<chainsN;j++)
         meansFile<<samplers[j]->getAverage(i).FF<<" ";
      meansFile<<sum2/chainsN<<endl;
   }
   meansFile.close();
   //}}}
   // Write thetaAct: {{{
   if(args.isSet("thetaActFileName")){
      ofstream actFile(args.getS("thetaActFileName").c_str());
      if(actFile.is_open()){
         actFile<<"# N "<<chainsN*min(samplesSave,samplesN)<<endl;
         for(j=0;j<chainsN;j++){
            for(i=0;i<Sof(samplers[j]->getThetaActLog());i++)
               actFile<<samplers[j]->getThetaActLog()[i]<<" ";
         }
         actFile<<endl;
         actFile.close();
      }else{
         message("WARNING: Main: Unable to write thetaAct log: %s.\n",(args.getS("thetaActFileName")).c_str());
      }
   }
   // }}}
   // Free memory: {{{
   for(j=0;j<chainsN;j++){
      delete samplers[j];
   }
//   delete [] samplers;
   //}}}
   message("Total samples: %ld\n",totalSamples);
}//}}}

extern "C" int estimateExpression(int *argc, char* argv[]) {//{{{
clearDataEE();
string programDescription =
"Estimates expression given precomputed probabilities of (observed) reads' alignments.\n\
   Uses MCMC sampling algorithm to produce relative abundance or RPKM.\n";
   buildTime(argv[0]);
   // Set options {{{
   ArgumentParser args;
   args.init(programDescription,"[prob file]",1);
   args.addOptionS("o","outFile","outFilePrefix",1,"Prefix for the output files.");
   args.addOptionS("O","outputType","outputType",0,"Output type (theta, RPKM, counts, tau).","counts");
   args.addOptionB("G","gibbs","gibbs",0,"Use gibbs sampling instead of collapsed gibbs sampling.");
   args.addOptionS("p","parFile","parFileName",0,"File containing parameters for the sampler, which can be otherwise specified by --MCMC* options. As the file is checked after every MCMC iteration, the parameters can be adjusted while running.");
   args.addOptionS("t","trInfoFile","trInfoFileName",0,"File containing transcript information. (Necessary for RPKM)");
   args.addOptionL("P","procN","procN",0,"Maximum number of threads to be used.",4);
   args.addOptionS("","thetaActFile","thetaActFileName",0,"File for logging noise parameter \theta^{act}.");
   args.addOptionL("","MCMC_burnIn","MCMC_burnIn",0,"Length of sampler's burn in period.",1000);
   args.addOptionL("","MCMC_samplesN","MCMC_samplesN",0,"Initial number of samples produced. Doubles after every iteration.",1000);
   args.addOptionL("","MCMC_samplesSave","MCMC_samplesSave",0,"Number of samples recorder for each chain at the end.",500);
   args.addOptionL("","MCMC_samplesNmax","MCMC_samplesNmax",0,"Maximum number of samples produced in one iteration. After producing samplesNmax samples sampler finishes.",50000);
   args.addOptionL("","MCMC_chainsN","MCMC_chainsN",0,"Number of parallel chains used. At least two chains will be used.",4);
   args.addOptionD("","MCMC_scaleReduction","MCMC_scaleReduction",0,"Target scale reduction, sampler finishes after this value is met.",1.2);
   args.addOptionD("","MCMC_dirAlpha","MCMC_dirAlpha",0,"Alpha parameter for the Dirichlet distribution.",1.0);
   if(!args.parse(*argc,argv))return 0;
   // }}}
   MyTimer timer;
#ifdef SUPPORT_OPENMP
   omp_set_num_threads(args.getL("procN"));
#endif

   gibbsParameters gPar;
//{{{ Initialization:

   gPar.setParameters(args);
   if(args.isSet("parFileName")){
      gPar.setParameters(args.getS("parFileName"));
   }
   
   if((args.getS("outputType") == "theta")||(args.getS("outputType") == "THETA")){
      outTypeI=THETA;
      outTypeS="theta";
   }else if((args.getS("outputType") == "RPKM")||(args.getS("outputType") == "rpkm")){
      outTypeI=RPKM;
      outTypeS="rpkm";
   }else if(args.getS("outputType") == "tau"){
      outTypeI=TAU;
      outTypeS="tau";
   }else{
      outTypeI=COVERAGE;
      outTypeS="counts";
      if(args.getS("outputType") != "counts")
         message("Using output type \"counts\"\n");
   }

   if(args.verbose)gPar.getAllParameters();

   //}}}
   // {{{ Read transcriptInfo and .prob file 
   if(args.verbose)message("Reading data.\n");
   if((!args.isSet("trInfoFileName"))||(!trInfo.readInfo(args.getS("trInfoFileName")))){
      if(outTypeI==RPKM){
         error("Main: Missing transcript info file. This will cause problems if producing RPKM.");
         return 1;
      }
      M=0;
   }else{
      M = trInfo.getM()+1;
   }
   readData(args);
   if(M<=0){
      error("Main: Invalid number of transcripts in .prob file.\n");
   }
   // }}}

   if(args.verbose)timer.split();
   if(args.verbose)message("Starting the sampler.\n");
   MCMC(gPar,args);
   // {{{ Transpose and merge sample file 
   if(transposeFiles(samplesFileNames,args.getS("outFilePrefix")+"."+outTypeS,args.verbose,failedMessage)){
      if(args.verbose)message("Sample files transposed. Deleting.\n");
      for(long i=0;i<Sof(samplesFileNames);i++){
         remove(samplesFileNames[i].c_str());
      }
   }else{
      message("Transposing files failed. Please check the files and try using trasposeLargeFile to transpose & merge the files into single file.\n");
   }
   //}}}
   if(args.verbose){message("DONE. "); timer.split(0,'h');}
   return 0;
}//}}}

#ifndef BIOC_BUILD
int main(int argc, char* argv[]) {
   return estimateExpression(&argc,argv);
}
#endif
