#include <R_ext/Rdynload.h>
#include <Rinternals.h>

extern "C" {
   int parseAlignment(int *argc,char* argv[]);
   int estimateExpression(int *argc, char* argv[]);
   int getVariance(int *argc, char* argv[]);
   int estimateHyperPar(int *argc, char* argv[]);
   int estimateDE(int *argc, char* argv[]);
}

static const R_CMethodDef cMethods[] = {
   {"_parseAlignment", (DL_FUNC) &parseAlignment, 2, (R_NativePrimitiveArgType[4]){INTSXP, STRSXP} },
   {"_estimateExpression", (DL_FUNC) &estimateExpression, 2, (R_NativePrimitiveArgType[4]){INTSXP, STRSXP} },
   {"_getVariance", (DL_FUNC) &getVariance, 2, (R_NativePrimitiveArgType[4]){INTSXP, STRSXP} },
   {"_estimateHyperPar", (DL_FUNC) &estimateHyperPar, 2, (R_NativePrimitiveArgType[4]){INTSXP, STRSXP} },
   {"_estimateDE", (DL_FUNC) &estimateDE, 2, (R_NativePrimitiveArgType[4]){INTSXP, STRSXP} },
   {NULL, NULL, 0}
};

extern "C" void R_init_BitSeq(DllInfo *info){
   R_registerRoutines(info, cMethods, NULL, NULL, NULL);

}
