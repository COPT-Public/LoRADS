
#include "def_asdp_func_set.h"
#include "asdp_algo.h"

extern void ASDPInitFuncSet(asdp_func **pfunc, int nLpCols){
    asdp_func *func;
    ASDP_INIT(func, asdp_func, 1);
    if (nLpCols > 0){
        // function with lp cone
        func->BMCalGrad = BMCalGradLP;
        func->BMCalq12p12 = BMCalq12p12LP;
        func->InitConstrValAll = ASDPInitConstrValAllLP;
        func->InitConstrValSum = ASDPInitConstrValSumLP;
        func->LBFGSDirUseGrad = LBFGSDirectionUseGradLP;
        func->LBFGSDirection = LBFGSDirectionLP;
        func->admmUpdateVar = ASDPUpdateSDPLPVar;
        func->calObj = ASDPCalObjLP;
        func->copyRtoV = copyRtoVLP;
        func->setAsNegGrad = SetyAsNegGradLP;
        func->setlbfgsHisTwo = setlbfgsHisTwoLP;
        func->updateDimac = ASDPUpdateDimacErrorLP;
        func->BMupdateVar = BMupdateVarLP;
    }else{
        // function with sdp cone only
        func->BMCalGrad = BMCalGrad;
        func->BMCalq12p12 = BMCalq12p12;
        func->InitConstrValAll = ASDPInitConstrValAll;
        func->InitConstrValSum = ASDPInitConstrValSum;
        func->LBFGSDirUseGrad = LBFGSDirectionUseGrad;
        func->LBFGSDirection = LBFGSDirection;
        func->admmUpdateVar = ASDPUpdateSDPVar;
        func->calObj = ASDPCalObj;
        func->copyRtoV = copyRtoV;
        func->setAsNegGrad = SetyAsNegGrad;
        func->setlbfgsHisTwo = setlbfgsHisTwo;
        func->updateDimac = ASDPUpdateDimacError;
        func->BMupdateVar = BMupdateVar;
    }
    *pfunc = func;
}
