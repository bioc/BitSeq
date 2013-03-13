#include <R_ext/Rdynload.h>
#include <Rinternals.h>

extern "C" {
   int parseAlignment(int *argc,char* argv[]);
   int estimateExpression(int *argc, char* argv[]);
   int getVariance(int *argc, char* argv[]);
   int estimateHyperPar(int *argc, char* argv[]);
   int estimateDE(int *argc, char* argv[]);
   int getGeneExpression(int *argc, char* argv[]);
   int getWithinGeneExpression(int *argc, char* argv[]);
}

static R_NativePrimitiveArgType my_types[] = {INTSXP, STRSXP};

static const R_CMethodDef cMethods[] = {
   {"_parseAlignment", (DL_FUNC) &parseAlignment, 2, my_types },
   {"_estimateExpression", (DL_FUNC) &estimateExpression, 2, my_types },
   {"_getVariance", (DL_FUNC) &getVariance, 2, my_types },
   {"_estimateHyperPar", (DL_FUNC) &estimateHyperPar, 2, my_types },
   {"_estimateDE", (DL_FUNC) &estimateDE, 2, my_types },
   {"_getGeneExpression", (DL_FUNC) &getGeneExpression, 2, my_types },
   {"_getWithinGeneExpression", (DL_FUNC) &getWithinGeneExpression, 2, my_types },
   {NULL, NULL, 0}
};

extern "C" void R_init_BitSeq(DllInfo *info){
   R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
