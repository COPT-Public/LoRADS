


#include "lorads_solver.h"
#include "lorads_vec_opts.h"
#include "lorads_alg_common.h"
int MAX_ALM_SUB_ITER;

extern void ALMSetGrad(lorads_solver *ASolver, lorads_sdp_cone *ACone, lorads_sdp_dense *R, lorads_sdp_dense *Grad, lorads_int iCone, double rho)
{
    LORADS_ZERO(Grad->matElem, double, Grad->nRows * Grad->rank);
    double *M1 = ASolver->var->M1temp;
    lorads_int n = ASolver->nRows;
    lorads_int incx = 1;
    if (iCone == 0)
    {
        LORADS_ZERO(M1, double, n);
        // M1 = -lambda + rho * (constrVal - b)
        /*for (lorads_int iConstr = 0; iConstr < ASolver->nRows; ++iConstr){
            M1[iConstr] = -ASolver->var->dualVar[iConstr] + ASolver->rho * (ASolver->rowRHS[iConstr] - ASolver->var->constrValSum[iConstr]);
        }*/
        double minusOne = -1.0;
        axpy(&n, &minusOne, ASolver->var->dualVar, &incx, M1, &incx);
        double negRho = -1 * rho;
        axpy(&n, &negRho, ASolver->rowRHS, &incx, M1, &incx);
        axpy(&n, &rho, ASolver->var->constrValSum, &incx, M1, &incx);
    }

    ACone->sdp_obj_sum->zeros(ACone->sdp_obj_sum->dataMat);
    ACone->addObjCoeff(ACone->coneData, ACone->sdp_obj_sum);
    ACone->sdpDataWSum(ACone->coneData, M1, ACone->sdp_obj_sum);

    // Note: obj mul_rk no row sparse method
    ACone->sdp_obj_sum->mul_rk(ACone->sdp_obj_sum->dataMat, R, Grad->matElem);
    lorads_int len = Grad->nRows * Grad->rank;
    double TWO = 2.0;
    scal(&len, &TWO, Grad->matElem, &incx);
}


extern void ALMCalGrad(lorads_solver *ASolver, lorads_lp_dense *rLpDummy, lorads_lp_dense *gradLpDummy, lorads_sdp_dense **R, lorads_sdp_dense **Grad, double *lagNormSquare, double rho)
{
    lorads_int incx = 1;
    double lagNorm2 = 0.0;
    lagNormSquare[0] = 0.0;
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        // S receive gradient info
        ALMSetGrad(ASolver, ASolver->SDPCones[iCone], R[iCone], Grad[iCone], iCone, rho);
        lorads_int n = Grad[iCone]->nRows * Grad[iCone]->rank;
        lagNorm2 = nrm2(&n, Grad[iCone]->matElem, &incx);
        lagNormSquare[0] += lagNorm2 * lagNorm2;
    }
}

extern void ALMSetGradLP(lorads_solver *ASolver, lorads_lp_cone *lp_cone, lorads_lp_dense *rLp, lorads_lp_dense *Grad, lorads_int iCol, double rho)
{
    Grad->matElem[iCol] = 0.0;
    double *M1 = ASolver->var->M1temp;
    lorads_int n = ASolver->nRows;
    lorads_int incx = 1;

    if (iCol == 0)
    {
        LORADS_ZERO(M1, double, n);
        // M1 = -laalmd - rho * b + rho * (sum A(uv))
        double minusOne = -1.0;
        axpy(&n, &minusOne, ASolver->var->dualVar, &incx, M1, &incx);
        double negRho = -1 * rho;
        axpy(&n, &negRho, ASolver->rowRHS, &incx, M1, &incx);
        axpy(&n, &rho, ASolver->var->constrValSum, &incx, M1, &incx);
    }

    double lpObjWlpDataSum = 0.0;
    lp_cone->objCoeffSum(lp_cone->coneData, &lpObjWlpDataSum, iCol);
    lp_cone->lpDataWSum(lp_cone->coneData, M1, &lpObjWlpDataSum, iCol);
    Grad->matElem[iCol] = 2 * lpObjWlpDataSum * rLp->matElem[iCol];
}


extern void ALMCalGradLP(lorads_solver *ASolver, lorads_lp_dense *rLp, lorads_lp_dense *gradLp, lorads_sdp_dense **R, lorads_sdp_dense **Grad, double *lagNormSquare, double rho)
{
    lorads_int incx = 1;
    double lagNorm2 = 0.0;
    lagNormSquare[0] = 0.0;
    for (lorads_int iCol = 0; iCol < ASolver->nLpCols; ++iCol)
    {
        ALMSetGradLP(ASolver, ASolver->lpCone, rLp, gradLp, iCol, rho);
    }
    lagNorm2 = nrm2(&(gradLp->nCols), gradLp->matElem, &incx);
    lagNormSquare[0] += lagNorm2 * lagNorm2;
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        // S receive gradient info
        ALMSetGrad(ASolver, ASolver->SDPCones[iCone], R[iCone], Grad[iCone], iCone, rho);
        lorads_int n = Grad[iCone]->nRows * Grad[iCone]->rank;
        lagNorm2 = nrm2(&n, Grad[iCone]->matElem, &incx);
        lagNormSquare[0] += lagNorm2 * lagNorm2;
    }
}

double LORADSnthroot(double base, lorads_int n) {
    if (base < 0 && n % 2 == 0) {
        return NAN;
    }
    if (base > 0){
        return pow(base, 1.0 / n);
    }else{
        return -pow(-base, 1.0 / n);
    }

}

extern lorads_int LORADScubic_equation(double a, double b, double c, double d, double *res){
    double A = b*b - 3*a*c;
    double B = b*c - 9*a*d;
    double C = c*c - 3*b*d;
    double delta = B*B - 4*A*C;
    res[0] = 0.0;
    res[1] = 0.0;
    res[2] = 0.0;
    if(A == 0 && B == 0){
        // one root
        res[0] = LORADS_MAX(res[0], -c/b);
        return 1;
    }else if (delta > 0){
        double Y1 = A*b + 1.5*a*(-B + sqrt(delta));
        double Y2 = A*b + 1.5*a*(-B - sqrt(delta));
        // double Y1_3 = Y1/(2*A*sqrt(A));
        double Y1_3 = LORADSnthroot(Y1, 3);
        double Y2_3 = LORADSnthroot(Y2, 3);
        res[0] = LORADS_MAX(res[0], (-b - Y1_3 - Y2_3)/3/a);
        return 1;
    }
    else if (delta == 0 && A != 0 && B != 0){
        double K = B/A;
        res[0] = -b/a + K;
        res[1] = -K/2;
        return 2;
    }else if (delta < 0){
        double sqA = sqrt(A);
        double T = (A*b - 1.5*a*B)/(A*sqA);
        double theta = acos(T);
        double csth = cos(theta/3);
        double sn3th = sqrt(3) * sin(theta/3);
        double root1 = (-b - 2*sqA*csth)/3/a;
        double root2 = (-b + sqA*(csth + sn3th))/3/a;
        double root3 = (-b + sqA*(csth - sn3th))/3/a;
        res[0] = root1; res[1] = root2; res[2] = root3;
        return 3;
    }else{
        return 0;
    }
}

static double functionVal(double a, double b, double c, double d, double x)
{
    return a * pow(x, 4) + b * pow(x, 3) + c * pow(x, 2) + d * x;
}

extern lorads_int ALMLineSearch(double rho, lorads_int n, double *lambd, double p1, double p2, double *q0, double *q1, double *q2, double *tau)
{
    lorads_int incx = 1;
    double q2nrm2 = nrm2(&n, q2, &incx);
    double a = rho * q2nrm2 * q2nrm2 / 2;
    double b = rho * dot(&n, q1, &incx, q2, &incx);
    // q0 = (1 / rho * lambd + q0)
    double rhoInv = 1 / rho;
    axpy(&n, &rhoInv, lambd, &incx, q0, &incx);
    double q1nrm = nrm2(&n, q1, &incx);
    double c = p2 - rho * dot(&n, q0, &incx, q2, &incx) + rho * q1nrm * q1nrm / 2;
    double d = p1 - rho * dot(&n, q0, &incx, q1, &incx);
    double roots[3] = {0.0, 0.0, 0.0};
    lorads_int rootNum = LORADScubic_equation(4 * a, 3 * b, 2 * c, d, roots);
    double f0 = 0.0;
    double tauMax = 1.0;
    double f1 = functionVal(a, b, c, d, 1.0);
    double fRoot1 = 1e+30;
    double fRoot2 = 1e+30;
    double fRoot3 = 1e+30;
    if (rootNum >= 1)
    {
        if (roots[0] > 1e-20 && roots[0] <= tauMax)
        {
            fRoot1 = functionVal(a, b, c, d, roots[0]);
        }
    }
    if (rootNum >= 2)
    {
        if (roots[1] > 1e-20 && roots[1] <= tauMax)
        {
            fRoot2 = functionVal(a, b, c, d, roots[1]);
        }
    }
    if (rootNum == 3)
    {
        if (roots[2] > 1e-20 && roots[2] <= tauMax)
        {
            fRoot3 = functionVal(a, b, c, d, roots[2]);
        }
    }
    else if (rootNum == 0)
    {
        ;
    }
    double minFval = LORADS_MIN(LORADS_MIN(LORADS_MIN(LORADS_MIN(f0, f1), fRoot1), fRoot2), fRoot3);
    if (fabs(minFval - f0) < 1e-10)
    {
        tau[0] = 0.0;
    }
    if (fabs(minFval - f1) < 1e-10)
    {
        tau[0] = 1.0;
    }
    if (fabs(minFval - fRoot1) < 1e-10)
    {
        tau[0] = roots[0];
    }
    if (fabs(minFval - fRoot2) < 1e-10)
    {
        tau[0] = roots[1];
    }
    if (fabs(minFval - fRoot3) < 1e-10)
    {
        tau[0] = roots[2];
    }
    return rootNum;
}

extern void LBFGSDirection(lorads_params *params, lorads_solver *ASolver, lbfgs_node *head, lorads_lp_dense *gradLpDummy, lorads_lp_dense *dDummy, lorads_sdp_dense **Grad, lorads_sdp_dense **D, lorads_int innerIter)
{
    /*head is oldest info, update head node and move head pointer lorads_into next*/
    lorads_int incx = 1;
    lorads_int n = 0;
    double minusOne = -1.0;

    lorads_int idx = 0;
    if (innerIter == 0){
        // D = -grad
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            LORADS_ZERO(D[iCone]->matElem, double, n);
            axpy(&n, &minusOne, Grad[iCone]->matElem, &incx, D[iCone]->matElem, &incx);
        }
        return;
    }else{
#ifdef FIX_INI_POINT
        double *Dtemp = ASolver->var->Dtemp;
        LORADS_ZERO(ASolver->var->Dtemp, double, head->allElem);
        idx = 0;
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            LORADS_MEMCPY(&Dtemp[idx], Grad[iCone]->matElem, double, n);
            idx += n;
//            if (iCone == 0){
//                printf("----------------\n");
//                for (int iElem = 0; iElem < 10; ++iElem){
//                    printf("Grad: (-D)[%d]: %f\n", iElem, Dtemp[iElem]);
//                }
//            }
        }

//        scal(&(head->allElem), &minusOne, Dtemp, &incx);
        lorads_int NodeNum = 0;
        if (innerIter <= params->lbfgsListLength - 1){
            NodeNum = innerIter;
        }else{
            NodeNum = params->lbfgsListLength;
        }
        lbfgs_node *node = head->prev;
        printf("NodeNum:%d", NodeNum);
        for (lorads_int k = 0; k < NodeNum; ++k){
            // new to old, head -> previous is newst
            // calculate alpha
            double temp = 0.0;
            temp = dot(&(head->allElem), node->s, &incx, Dtemp, &incx);
//            double nrm_s = nrm2(&(head->allElem), node->s, &incx);
//            double nrm_D = nrm2(&(head->allElem), Dtemp, &incx);
//            printf("nrm_s: %f\n", nrm_s);
//            printf("nrm_D: %f\n", nrm_D);
//            lorads_int temp_n = 10000;
//            double temp10 = dot(&temp_n, node->s, &incx, Dtemp, &incx);
//            for (int iElem = 0; iElem < 5; iElem ++){
//                printf("**node->s[%d]:%f\n", iElem, node->s[iElem]);
//                printf("**Dtemp[%d]:%f\n", iElem, Dtemp[iElem]);
//            }
//            printf("beta: %.20f\n", node->beta);
//            printf("temp: %f\n", temp);
            node->alpha = node->beta * temp;
//            printf("alpha: %f\n", node->alpha);
            double alphaNeg = -1 * node->alpha;
            axpy(&(head->allElem), &alphaNeg, node->y, &incx, Dtemp, &incx);
            node = node->prev;
        }
//        printf("----------------\n");
        node = node->next;
        for (lorads_int k = 0; k < NodeNum; ++k){
            double weight = node->alpha - node->beta * dot(&(head->allElem), node->y, &incx, Dtemp, &incx);
            axpy(&(head->allElem), &weight, node->s, &incx, Dtemp, &incx);
            node = node->next;
        }
        scal(&(head->allElem), &minusOne, Dtemp, &incx);
        idx = 0;
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            LORADS_MEMCPY(D[iCone]->matElem, &Dtemp[idx], double, n);
            idx += n;
        }
#else
//        double *Dtemp = ASolver->var->Dtemp;
//        LORADS_ZERO(ASolver->var->Dtemp, double, head->allElem);
//        idx = 0;
//        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
//            n = Grad[iCone]->nRows * Grad[iCone]->rank;
//            LORADS_MEMCPY(&Dtemp[idx], Grad[iCone]->matElem, double, n);
//            idx += n;
//        }
//        lorads_int NodeNum = 0;
//        if (innerIter <= params->lbfgsListLength){
//            NodeNum = innerIter;
//        }else{
//            NodeNum = params->lbfgsListLength;
//        }
//        lbfgs_node *node = head->prev;
//        for (lorads_int k = 0; k < NodeNum; ++k){
//            // new to old, head -> previous is newst
//            // calculate alpha
//            double temp = 0.0;
//            temp = dot(&(head->allElem), node->s, &incx, Dtemp, &incx);
//
//            node->alpha = node->beta * temp;
//            double alphaNeg = -1 * node->alpha;
//            axpy(&(head->allElem), &alphaNeg, node->y, &incx, Dtemp, &incx);
//            node = node->prev;
//        }
//        node = node->next;
//        for (lorads_int k = 0; k < NodeNum; ++k){
//            double weight = node->alpha - node->beta * dot(&(head->allElem), node->y, &incx, Dtemp, &incx);
//            axpy(&(head->allElem), &weight, node->s, &incx, Dtemp, &incx);
//            node = node->next;
//        }
//        scal(&(head->allElem), &minusOne, Dtemp, &incx);
//        idx = 0;
//        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
//            n = Grad[iCone]->nRows * Grad[iCone]->rank;
//            LORADS_MEMCPY(D[iCone]->matElem, &Dtemp[idx], double, n);
//            idx += n;
//        }
//        return;
        double *Dtemp = ASolver->var->Dtemp;
        LORADS_ZERO(ASolver->var->Dtemp, double, head->allElem);
        idx = 0;
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            LORADS_MEMCPY(&Dtemp[idx], Grad[iCone]->matElem, double, n);
            idx += n;
        }
        lorads_int NodeNum = 0;
        if (innerIter <= params->lbfgsListLength - 1){
            NodeNum = innerIter;
        }else{
            NodeNum = params->lbfgsListLength;
        }
        lbfgs_node *node = head->prev;
        for (lorads_int k = 0; k < NodeNum; ++k){
            // new to old, head -> previous is newst
            // calculate alpha
            double temp = 0.0;
            temp = dot(&(head->allElem), node->s, &incx, Dtemp, &incx);
            node->alpha = node->beta * temp;
            double alphaNeg = -1 * node->alpha;
            axpy(&(head->allElem), &alphaNeg, node->y, &incx, Dtemp, &incx);
            node = node->prev;
        }
        node = node->next;
        for (lorads_int k = 0; k < NodeNum; ++k){
            double weight = node->alpha - node->beta * dot(&(head->allElem), node->y, &incx, Dtemp, &incx);
            axpy(&(head->allElem), &weight, node->s, &incx, Dtemp, &incx);
            node = node->next;
        }
        scal(&(head->allElem), &minusOne, Dtemp, &incx);
        idx = 0;
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            LORADS_MEMCPY(D[iCone]->matElem, &Dtemp[idx], double, n);
            idx += n;
        }
#endif
    }
}

extern void LBFGSDirectionLP(lorads_params *params, lorads_solver *ASolver, lbfgs_node *head, lorads_lp_dense *gradLp, lorads_lp_dense *d, lorads_sdp_dense **Grad, lorads_sdp_dense **D, lorads_int innerIter)
{
    /*head is oldest info, update head node and move head pointer lorads_into next*/
    lorads_int incx = 1;
    lorads_int n = 0;
    double minusOne = -1.0;
    double *Dtemp;
    lorads_int idx = 0;
    if (innerIter == 0)
    {
        // D = -grad
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            LORADS_ZERO(D[iCone]->matElem, double, n);
            axpy(&n, &minusOne, Grad[iCone]->matElem, &incx, D[iCone]->matElem, &incx);
        }
        LORADS_ZERO(d->matElem, double, d->nCols);
        n = d->nCols;
        axpy(&n, &minusOne, gradLp->matElem, &incx, d->matElem, &incx);
        return;
    }
    else
    {
        // Dtemp = grad
        LORADS_INIT(Dtemp, double, head->allElem);
        idx = 0;
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            LORADS_MEMCPY(&Dtemp[idx], Grad[iCone]->matElem, double, n);
            idx += n;
        }
        LORADS_MEMCPY(&Dtemp[idx], gradLp->matElem, double, gradLp->nCols);
    }

    lorads_int NodeNum = 0;
    if (innerIter <= params->lbfgsListLength)
    {
        NodeNum = innerIter;
    }
    else
    {
        NodeNum = params->lbfgsListLength;
    }
    lbfgs_node *node = head->prev;
    for (lorads_int k = 0; k < NodeNum; ++k)
    {
        double temp = 0.0;
        temp = dot(&(head->allElem), node->s, &incx, Dtemp, &incx);
        node->alpha = node->beta * temp;
        double alphaNeg = -1 * node->alpha;
        axpy(&(head->allElem), &alphaNeg, node->y, &incx, Dtemp, &incx);
        node = node->prev;
    }
    node = node->next;
    for (lorads_int k = 0; k < NodeNum; ++k)
    {
        double weight = node->alpha - node->beta * dot(&(head->allElem), node->y, &incx, Dtemp, &incx);
        axpy(&(head->allElem), &weight, node->s, &incx, Dtemp, &incx);
        node = node->next;
    }
    scal(&(head->allElem), &minusOne, Dtemp, &incx);
    idx = 0;
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        n = Grad[iCone]->nRows * Grad[iCone]->rank;
        LORADS_MEMCPY(D[iCone]->matElem, &Dtemp[idx], double, n);
        idx += n;
    }
    LORADS_MEMCPY(d->matElem, &Dtemp[idx], double, d->nCols);
    idx += d->nCols;
    LORADS_FREE(Dtemp);
    return;
}

extern void LBFGSDirectionUseGrad(lorads_solver *ASolver, lorads_lp_dense *dDummy, lorads_lp_dense *gradLPDummy, lorads_sdp_dense **D, lorads_sdp_dense **Grad)
{
    double innerProduct = 0.0;
    lorads_int nElem = 0;
    lorads_int incx = 1;
    double minusOne = -1.0;
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        nElem = D[iCone]->nRows * D[iCone]->rank;
        innerProduct += dot(&nElem, D[iCone]->matElem, &incx, Grad[iCone]->matElem, &incx);
    }
    if (innerProduct >= 0)
    {
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            nElem = D[iCone]->nRows * D[iCone]->rank;
            LORADS_MEMCPY(D[iCone]->matElem, Grad[iCone]->matElem, double, nElem);
            scal(&nElem, &minusOne, D[iCone]->matElem, &incx);
        }
    }
}

extern void LBFGSDirectionUseGradLP(lorads_solver *ASolver, lorads_lp_dense *d, lorads_lp_dense *gradLP, lorads_sdp_dense **D, lorads_sdp_dense **Grad)
{
    double innerProduct = 0.0;
    lorads_int nElem = 0;
    lorads_int incx = 1;
    double minusOne = -1.0;
    innerProduct += dot(&(gradLP->nCols), d->matElem, &incx, gradLP->matElem, &incx);
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        nElem = D[iCone]->nRows * D[iCone]->rank;
        innerProduct += dot(&nElem, D[iCone]->matElem, &incx, Grad[iCone]->matElem, &incx);
    }
    if (innerProduct >= 0)
    {
        LORADS_MEMCPY(d->matElem, gradLP->matElem, double, gradLP->nCols);
        scal(&(gradLP->nCols), &minusOne, d->matElem, &incx);
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            nElem = D[iCone]->nRows * D[iCone]->rank;
            LORADS_MEMCPY(D[iCone]->matElem, Grad[iCone]->matElem, double, nElem);
            scal(&nElem, &minusOne, D[iCone]->matElem, &incx);
        }
    }
}

extern void LORADSConstrValSumALMtemp(lorads_solver *ASolver, double *q1)
{
    double alpha = 1.0;
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASolver->var->constrVal[iCone]->add(&alpha, ASolver->var->constrVal[iCone]->data, q1);
    }
}

extern void LORADSConstrValSumALMtempLP(lorads_solver *ASolver, double *q1)
{
    double alpha = 1.0;
    for (lorads_int iCol = 0; iCol < ASolver->nLpCols; ++iCol)
    {
        ASolver->var->constrValLP[iCol]->add(&alpha, ASolver->var->constrValLP[iCol]->data, q1);
    }
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASolver->var->constrVal[iCone]->add(&alpha, ASolver->var->constrVal[iCone]->data, q1);
    }
}



extern void ALMCalq12p12(lorads_solver *ASolver, lorads_lp_dense *rLpDummy, lorads_lp_dense *dDummy, lorads_sdp_dense **R, lorads_sdp_dense **D, double *q1, double *q2, double *p12)
{
    lorads_int n = ASolver->nRows;
    double one = 1.0;
    double zero = 0.0;
    lorads_int incx = 1;
    // LORADS_ZERO(q1, double, n);
    // LORADS_ZERO(q2, double, n);
    scal(&n, &zero, q1, &incx);
    scal(&n, &zero, q2, &incx);
    p12[0] = 0.0;
    p12[1] = 0.0;
    LORADSObjConstrValAll(ASolver, R, D, &p12[0]);
    LORADSConstrValSumALMtemp(ASolver, q1);
    double Two = 2.0;
    scal(&n, &Two, q1, &incx);
    p12[0] *= 2;

    LORADSObjConstrValAll(ASolver, D, D, &p12[1]);
    LORADSConstrValSumALMtemp(ASolver, q2);
}


extern void ALMCalq12p12LP(lorads_solver *ASolver, lorads_lp_dense *rLp, lorads_lp_dense *d, lorads_sdp_dense **R, lorads_sdp_dense **D, double *q1, double *q2, double *p12)
{
    lorads_int n = ASolver->nRows;
    double one = 1.0;
    lorads_int incx = 1;
    LORADS_ZERO(q1, double, n);
    LORADS_ZERO(q2, double, n);
    p12[0] = 0.0;
    p12[1] = 0.0;
    LORADSObjConstrValAllLP(ASolver, rLp, d, R, D, &p12[0]);
    LORADSConstrValSumALMtempLP(ASolver, q1);
    double Two = 2.0;
    scal(&n, &Two, q1, &incx);
    p12[0] *= 2;

    LORADSObjConstrValAllLP(ASolver, d, d, D, D, &p12[1]);
    LORADSConstrValSumALMtempLP(ASolver, q2);
}


extern void SetyAsNegGrad(lorads_solver *ASolver, lorads_lp_dense *gradLpDummy, lorads_sdp_dense **Grad)
{
    double minusOne = -1.0;
    lorads_int incx = 1;
    lorads_int n = 0;
    lorads_int idx = 0;
    lbfgs_node *head = ASolver->lbfgsHis;
    LORADS_ZERO(head->y, double, head->allElem);
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        lorads_sdp_dense *gradICone = Grad[iCone];
        n = gradICone->nRows * gradICone->rank;
        axpy(&n, &minusOne, gradICone->matElem, &incx, &head->y[idx], &incx);
        idx += n;
    }
}

extern void SetyAsNegGradLP(lorads_solver *ASolver, lorads_lp_dense *gradLp, lorads_sdp_dense **Grad)
{
    double minusOne = -1.0;
    lorads_int incx = 1;
    lorads_int n = 0;
    lorads_int idx = 0;
    lbfgs_node *head = ASolver->lbfgsHis;
    LORADS_ZERO(head->y, double, head->allElem);
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        lorads_sdp_dense *gradICone = Grad[iCone];
        n = gradICone->nRows * gradICone->rank;
        axpy(&n, &minusOne, gradICone->matElem, &incx, &head->y[idx], &incx);
        idx += n;
    }
    axpy(&gradLp->nCols, &minusOne, gradLp->matElem, &incx, &head->y[idx], &incx);
    // idx += gradLp->nCols;
}

extern void ALMupdateVar(lorads_solver *ASolver, lorads_lp_dense *rLpDummy, lorads_lp_dense *dDummy, lorads_sdp_dense **R, lorads_sdp_dense **D, double tau)
{
    lorads_int incx = 1;
    lorads_int nElem = 0;
#ifdef FIX_INI_POINT
    double check = 0.0;
    double temp = 0.0;
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        nElem = R[iCone]->nRows * R[iCone]->rank;
        axpy(&nElem, &tau, D[iCone]->matElem, &incx, R[iCone]->matElem, &incx);
//        if (iCone == 0){
//            for (int iElem = 0; iElem < 5; ++iElem){
//                printf("tau * D[%d]: %f\n", iElem, tau * D[iCone]->matElem[iElem]);
//        }
//        }

        temp = nrm2(&nElem, R[iCone]->matElem, &incx);
        check += temp * temp;
    }
        printf("check: %f\n", sqrt(check));
#else
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        nElem = R[iCone]->nRows * R[iCone]->rank;
        axpy(&nElem, &tau, D[iCone]->matElem, &incx, R[iCone]->matElem, &incx);
    }
#endif

}

extern void ALMupdateVarLP(lorads_solver *ASolver, lorads_lp_dense *rLp, lorads_lp_dense *d, lorads_sdp_dense **R, lorads_sdp_dense **D, double tau)
{
    lorads_int incx = 1;
    ALMupdateVar(ASolver, rLp, d, R, D, tau);
    axpy(&(rLp->nCols), &tau, d->matElem, &incx, rLp->matElem, &incx);
}

extern void setlbfgsHisTwo(lorads_solver *ASolver, lorads_lp_dense *gradLpDummy, lorads_lp_dense *dDummy, lorads_sdp_dense **Grad, lorads_sdp_dense **D, double tau)
{
    double minusOne = -1.0;
    double One = 1.0;
    lorads_int incx = 1;
    lorads_int n = 0;
    lorads_int idx = 0;
    lbfgs_node *head = ASolver->lbfgsHis;
    LORADS_ZERO(head->s, double, head->allElem);
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        lorads_sdp_dense *gradICone = Grad[iCone];
        n = gradICone->nRows * gradICone->rank;
        // head->s = tau * D
        axpy(&n, &tau, D[iCone]->matElem, &incx, &head->s[idx], &incx);
        // head->y += GradNew
        axpy(&n, &One, Grad[iCone]->matElem, &incx, &head->y[idx], &incx);
        idx += n;
    }
    head->beta = 1.0 / dot(&(head->allElem), head->y, &incx, head->s, &incx);
    ASolver->lbfgsHis = head->next;
}

extern void setlbfgsHisTwoLP(lorads_solver *ASolver, lorads_lp_dense *gradLp, lorads_lp_dense *d, lorads_sdp_dense **Grad, lorads_sdp_dense **D, double tau)
{
    double minusOne = -1.0;
    double One = 1.0;
    lorads_int incx = 1;
    lorads_int n = 0;
    lorads_int idx = 0;
    lbfgs_node *head = ASolver->lbfgsHis;
    LORADS_ZERO(head->s, double, head->allElem);
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        lorads_sdp_dense *gradICone = Grad[iCone];
        n = gradICone->nRows * gradICone->rank;
        // head->s = tau * D
        axpy(&n, &tau, D[iCone]->matElem, &incx, &head->s[idx], &incx);
        // head->y += GradNew
        axpy(&n, &One, Grad[iCone]->matElem, &incx, &head->y[idx], &incx);
        idx += n;
    }
    axpy(&(d->nCols), &tau, d->matElem, &incx, &head->s[idx], &incx);
    axpy(&(d->nCols), &One, gradLp->matElem, &incx, &head->y[idx], &incx);
    // idx += d->nCols;
    head->beta = 1.0 / dot(&(head->allElem), head->y, &incx, head->s, &incx);
    ASolver->lbfgsHis = head->next;
}


static void rootNum0PrintInfo(lorads_int minIter){
#ifdef LORADS_INT32
    printf("*Numerical Fail in ALM minter %d, we use the best feasible solution as warm start\n", minIter);
#endif
#ifdef UNIX_INT64
    printf("*Numerical Fail in ALM minter %ld, we use the best feasible solution as warm start\n", minIter);
#endif
#ifdef MAC_INT64
    printf("*Numerical Fail in ALM minter %lld, we use the best feasible solution as warm start\n", minIter);
#endif
}

static void ALMPrintLog(lorads_alm_state *alm_state_pointer, double time){
#ifdef INT32
    printf("ALM OuterIter:%d InnerIter:%d pObj:%5.5e dObj:%5.5e pInfea(1):%5.5e pInfea(Inf):%5.5e pdGap:%5.5e rho:%3.2f Time:%3.2f\n",
           alm_state_pointer->outerIter, alm_state_pointer->innerIter,
           alm_state_pointer->primal_objective_value, alm_state_pointer->dual_objective_value,
           alm_state_pointer->l_1_primal_infeasibility, alm_state_pointer->l_inf_primal_infeasibility,
           alm_state_pointer->primal_dual_gap, alm_state_pointer->rho, time);
#endif

#ifdef UNIX_INT64
    printf("OuterIter:%ld InnerIter:%ld pObj:%5.5e dObj:%5.5e pInfea(1):%5.5e pInfea(Inf):%5.5e pdGap:%5.5e rho:%3.2f Time:%3.2f\n",
           alm_state_pointer->outerIter, alm_state_pointer->innerIter,
           alm_state_pointer->primal_objective_value, alm_state_pointer->dual_objective_value,
           alm_state_pointer->l_1_primal_infeasibility, alm_state_pointer->l_inf_primal_infeasibility,
           alm_state_pointer->primal_dual_gap, alm_state_pointer->rho, time);
#endif

#ifdef MAC_INT64
    printf("OuterIter:%lld InnerIter:%lld pObj:%5.5e dObj:%5.5e pInfea(1):%5.5e pInfea(Inf):%5.5e pdGap:%5.5e rho:%3.2f Time:%3.2f\n",
           alm_state_pointer->outerIter, alm_state_pointer->innerIter,
           alm_state_pointer->primal_objective_value, alm_state_pointer->dual_objective_value,
           alm_state_pointer->l_1_primal_infeasibility, alm_state_pointer->l_inf_primal_infeasibility,
           alm_state_pointer->primal_dual_gap, alm_state_pointer->rho, time);
#endif
}

extern lorads_int LORADS_ALMOptimize_reopt(lorads_params *params, lorads_solver *ASolver, lorads_alm_state *alm_iter_state, bool early_stop, double rho_update_factor, double timeSolveStart)
{
    double ori_start = LUtilGetTimeStamp();
    bool is_rank_max = CheckAllRankMax(ASolver, 1.0);
    lorads_func *aFunc;
    LORADSInitFuncSet(&aFunc, ASolver->nLpCols);
    lorads_int retcode = LORADS_RETCODE_OK;
    lorads_int incx = 1;
    lorads_int last_outter_iter_start = 1;
    double minusOne = -1.0;
    double tau = 0.0;
    double rho_certificate, rho_certificate_tol, rho_certificate_val, lagNormSquare, bestInfe;
    ALG_START:
    rho_certificate = 0.1;
    rho_certificate_tol = rho_certificate / alm_iter_state->rho;
    rho_certificate_val = 0;
    lagNormSquare = 0.0;
    bestInfe = 1e+30;
    aFunc->InitConstrValAll(ASolver, ASolver->var->rLp, ASolver->var->rLp, ASolver->var->R, ASolver->var->R);
    aFunc->InitConstrValSum(ASolver);

    aFunc->ALMCalGrad(ASolver, ASolver->var->rLp, ASolver->var->gradLp, ASolver->var->R, ASolver->var->Grad, &lagNormSquare, alm_iter_state->rho);
    rho_certificate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
    char difficulty = HARD;
    lorads_int localIter = 0;
    lorads_int clearLBFGS = 0;
    lorads_int rank_flag = 0;
    double rank_update_factor = 1.5;
    lorads_int k = alm_iter_state->outerIter;
    lorads_int k0 = alm_iter_state->outerIter;
    lorads_int rho_factor_flag = 0;
    lorads_int rank_flag_thres;
    if (params->dyrankLevel == 0){
        rank_flag_thres = 1e8;
    }else if (params->dyrankLevel == 1){
        rank_flag_thres = 150;
    }else if (params->dyrankLevel == 2){
        rank_flag_thres = 15;
    } else if (params->dyrankLevel == 3){
        rank_flag_thres = 5;
    }

    int max_sub_iter_inc_factor = 10000;
    int max_sub_iter_ceil = 25000;
    int update_max_sub_iter_counter = 0;

    while (true)
    {
        if ((k > params->maxALMIter) && ((alm_iter_state->l_inf_primal_infeasibility <= params->phase1Tol) && ((alm_iter_state->primal_dual_gap <= LORADS_MAX(params->phase1Tol, params->phase2Tol * 5)) || (!params->highAccMode))))
        {
            break;
        }

        lorads_int rho_inner_iter = 0;
        double rank_update_alpha = 0.1;
        double rank_update_threshold = 0.005;
        double rank_update_current_ema = 0.0;
        double rank_update_old_ema = 0.0;
        lorads_int rank_update_evaluate_interval = 5;
        lorads_int rank_update_counter = 1;
        lorads_int cur_iter_counter = 1;
        if (update_max_sub_iter_counter >= 2){
            update_max_sub_iter_counter = 0;
            MAX_ALM_SUB_ITER += max_sub_iter_inc_factor;
            MAX_ALM_SUB_ITER = LORADS_MIN(MAX_ALM_SUB_ITER, max_sub_iter_ceil);
        }
        while (difficulty != EASY)
        {
            localIter = 0;
            lorads_int if_break = LUtilUpdateCheckEma(&rank_update_current_ema, &rank_update_old_ema, rho_certificate_val, rank_update_alpha, rank_update_threshold, rank_update_evaluate_interval, &rank_update_counter);
            if (!if_break && !params->highAccMode)
            {
                break;
            }
            if (cur_iter_counter >= MAX_ALM_SUB_ITER)
            {
                update_max_sub_iter_counter += 1;
                break;
            }
            if (rank_flag >= rank_flag_thres && (!is_rank_max) && (k - last_outter_iter_start >= 3))
            {
                break;
            }
            if (rho_certificate_val <= rho_certificate_tol)
            {
                break;
            }
            while (rho_certificate_val - rho_certificate_tol > params->endALMSubTol)
            {
                if ((localIter - 1) % 300 == 0)
                {
                    clearLBFGS = 0;
                }
                aFunc->LBFGSDirection(params, ASolver, ASolver->lbfgsHis, ASolver->var->gradLp, ASolver->var->uLp, ASolver->var->Grad, ASolver->var->U, clearLBFGS);
                aFunc->LBFGSDirUseGrad(ASolver, ASolver->var->uLp, ASolver->var->gradLp, ASolver->var->U, ASolver->var->Grad);
                double *q0 = ASolver->var->M1temp;
                // q0 = b - A(RRt)
                LORADS_MEMCPY(q0, ASolver->rowRHS, double, ASolver->nRows);
                axpy(&(ASolver->nRows), &minusOne, ASolver->var->constrValSum, &incx, q0, &incx);

                double p12[2];
                aFunc->ALMCalq12p12(ASolver, ASolver->var->rLp, ASolver->var->uLp, ASolver->var->R, ASolver->var->U, ASolver->var->ARDSum, ASolver->var->ADDSum, p12);
                lorads_int rootNum = ALMLineSearch(alm_iter_state->rho, ASolver->nRows, ASolver->var->dualVar, p12[0], p12[1], q0, ASolver->var->ARDSum, ASolver->var->ADDSum, &tau);
                if (rootNum == 0){
                    retcode = RET_CODE_NUM_ERR;
                    goto END_ALM;
                }
                if (fabs(tau) < params->endTauTol){
                    printf("update rho, tau is too small :%5.3e\n", tau);
                    alm_iter_state->innerIter++;
                    localIter++;
                    cur_iter_counter++;
                    clearLBFGS++;
                    goto UpdateRho;
                }
                // y = (-grad)
                // set lbfgs one
                aFunc->setAsNegGrad(ASolver, ASolver->var->gradLp, ASolver->var->Grad);
                // update R
                aFunc->ALMupdateVar(ASolver, ASolver->var->rLp, ASolver->var->uLp, ASolver->var->R, ASolver->var->U, tau);
                // update gradient
                lagNormSquare = 0.0;
                // update constrValSum first
                double tauSquare = tau * tau;
                axpy(&(ASolver->nRows), &tau, ASolver->var->ARDSum, &incx, ASolver->var->constrValSum, &incx);
                axpy(&(ASolver->nRows), &(tauSquare), ASolver->var->ADDSum, &incx, ASolver->var->constrValSum, &incx);
                aFunc->ALMCalGrad(ASolver, ASolver->var->rLp, ASolver->var->gradLp, ASolver->var->R, ASolver->var->Grad, &lagNormSquare, alm_iter_state->rho);
                rho_certificate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
                aFunc->setlbfgsHisTwo(ASolver, ASolver->var->gradLp, ASolver->var->uLp, ASolver->var->Grad, ASolver->var->U, tau);

                aFunc->updateDimacsALM(ASolver, ASolver->var->R, ASolver->var->R, ASolver->var->rLp, ASolver->var->rLp);
                alm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
                alm_iter_state->l_inf_primal_infeasibility = alm_iter_state->l_1_primal_infeasibility * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf);
                alm_iter_state->innerIter++;
                localIter++;
                cur_iter_counter++;
                clearLBFGS++;
                if (localIter > 800){
                    break;
                }
            }
            LORADSUpdateDualVar(ASolver, alm_iter_state->rho);
            aFunc->ALMCalGrad(ASolver, ASolver->var->rLp, ASolver->var->gradLp, ASolver->var->R, ASolver->var->Grad, &lagNormSquare, alm_iter_state->rho);
            rho_certificate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
            if (localIter <= 20){
                difficulty = EASY;
            }
            else if (20 < localIter && localIter <= 100){
                difficulty = MEDIUM;
                rank_flag += 2;
            }
            else if (100 < localIter){
                difficulty = HARD;
                rank_flag += 3;
            }
            else if (400 < localIter){
                difficulty = SUPER;
                rank_flag += 4;
            }
            if (difficulty == EASY){
                rank_flag = 0;
            }
            rho_inner_iter ++;
        }
        UpdateRho:
        do{
            alm_iter_state->rho *= rho_update_factor;
            aFunc->ALMCalGrad(ASolver, ASolver->var->rLp, ASolver->var->gradLp, ASolver->var->R, ASolver->var->Grad, &lagNormSquare, alm_iter_state->rho);
            rho_certificate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
            rho_certificate_tol = rho_certificate / alm_iter_state->rho;
        }while(rho_certificate_tol >= rho_certificate_val);
        if (alm_iter_state->rho >= 5e4 && rho_factor_flag < 4){
            rho_update_factor = sqrt(sqrt(rho_update_factor));
            rho_factor_flag = 4;
        }else if(alm_iter_state->rho >= 5e6 && rho_factor_flag < 6){
            rho_update_factor = sqrt(sqrt(rho_update_factor));
            rho_factor_flag = 6;
        }else if (alm_iter_state->rho >= 5e8 && rho_factor_flag < 8){
            rho_update_factor = sqrt(sqrt(rho_update_factor));
            rho_factor_flag = 8;
        }
        difficulty = HARD;
        clearLBFGS = 0;
        k += 1;
        alm_iter_state->outerIter = k;
        if (k % 1 == 0){
            aFunc->calObj_alm(ASolver);
            LORADSCalDualObj(ASolver);
            aFunc->updateDimacsALM(ASolver, ASolver->var->R, ASolver->var->R, ASolver->var->rLp, ASolver->var->rLp);
            alm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
            alm_iter_state->primal_objective_value = ASolver->pObjVal;
            alm_iter_state->dual_objective_value = ASolver->dObjVal;
            alm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
            alm_iter_state->l_inf_primal_infeasibility = alm_iter_state->l_1_primal_infeasibility * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf);
            alm_iter_state->l_1_dual_infeasibility = 99; // does not evaluate dual error
            alm_iter_state->l_inf_dual_infeasibility = 99; // does not evaluate dual error
            if (early_stop)
            {
                if (alm_iter_state->l_1_primal_infeasibility <= params->phase1Tol && alm_iter_state->primal_dual_gap <= LORADS_MAX(params->phase1Tol, params->phase2Tol * 5) && (k - k0) > 1)
                {
                    goto PRINT_AND_EXIT;
                }
                // if ((k > params->maxALMIter) && alm_iter_state->l_inf_primal_infeasibility <= params->phase2Tol && (k - k0) > 1)
                // {
                //     goto PRINT_AND_EXIT;
                // }
            }else{
                if(alm_iter_state->primal_dual_gap <= params->phase2Tol && alm_iter_state->l_1_primal_infeasibility <= params->phase2Tol && (k - k0) > 1){
                    goto PRINT_AND_EXIT;
                }
            }
            ALMPrintLog(alm_iter_state, LUtilGetTimeStamp() - ori_start);
            if (LUtilGetTimeStamp() - timeSolveStart >= params->timeSecLimit){
                goto PRINT_AND_EXIT;
            }
        }
        if (rank_flag >= rank_flag_thres && !is_rank_max && (ASolver->nCones <= 10)){
            rank_flag = 0;
            if (k - last_outter_iter_start >= 2 && ASolver->nCones <= 1e+10)
            {
                printf("increase the rank, factor:%f.\n", rank_update_factor);
                is_rank_max = AUG_RANK(ASolver, ASolver->var->rankElem, ASolver->nCones, rank_update_factor);
                alm_iter_state->outerIter = k;
                last_outter_iter_start = alm_iter_state->outerIter;
                goto ALG_START;
            }
        }
    }
    END_ALM:
    aFunc->calObj_alm(ASolver);
    LORADSCalDualObj(ASolver);
    aFunc->updateDimacsALM(ASolver, ASolver->var->R, ASolver->var->R, ASolver->var->rLp, ASolver->var->rLp);
    alm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
    alm_iter_state->l_1_primal_infeasibility = alm_iter_state->l_inf_primal_infeasibility * (1 + ASolver->bRHSNrmInf) / (1 + ASolver->bRHSNrm1);
    alm_iter_state->l_1_dual_infeasibility = 99; // does not evaluate dual error
    alm_iter_state->l_inf_dual_infeasibility = 99; // does not evaluate dual error
    PRINT_AND_EXIT:
    printf("-----------------------------------------------------------------------\n");
    printf("Exit ALM:\n");
    ALMPrintLog(alm_iter_state, LUtilGetTimeStamp() - ori_start);
    printf("-----------------------------------------------------------------------\n");
    return retcode;
}



extern lorads_int LORADS_ALMOptimize(lorads_params *params, lorads_solver *ASolver, lorads_alm_state *alm_iter_state, double rho_update_factor, double timeSolveStart)
{
    MAX_ALM_SUB_ITER = 5000;
    double ori_start = LUtilGetTimeStamp();
    bool is_rank_max = CheckAllRankMax(ASolver, 1.0);
    lorads_func *aFunc;
    LORADSInitFuncSet(&aFunc, ASolver->nLpCols);
    lorads_int retcode = LORADS_RETCODE_OK;
    lorads_int incx = 1;
    lorads_int last_outter_iter_start = 1;
    double minusOne = -1.0;
    double tau = 0.0;
    double rho_certificate, rho_certificate_tol, rho_certificate_val, lagNormSquare, bestInfe;
    ALG_START:
    rho_certificate = 0.1;
    rho_certificate_tol = rho_certificate / alm_iter_state->rho;
    rho_certificate_val = 0;
    lagNormSquare = 0.0;
    bestInfe = 1e+30;
    aFunc->InitConstrValAll(ASolver, ASolver->var->rLp, ASolver->var->rLp, ASolver->var->R, ASolver->var->R);
    aFunc->InitConstrValSum(ASolver);

    aFunc->ALMCalGrad(ASolver, ASolver->var->rLp, ASolver->var->gradLp, ASolver->var->R, ASolver->var->Grad, &lagNormSquare, alm_iter_state->rho);
    rho_certificate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
    char difficulty = HARD;
    lorads_int localIter = 0;
    lorads_int clearLBFGS = 0;
    lorads_int rank_flag = 0;
    double rank_update_factor = 1.5;
    rho_update_factor = params->ALMRhoFactor;
    lorads_int rho_factor_flag = 0;
    double rank_flag_thres;
    if (params->dyrankLevel == 0){
        rank_flag_thres = 1e8;
    }else if(params->dyrankLevel == 1){
        rank_flag_thres = 150;
    }else if(params->dyrankLevel == 2){
        rank_flag_thres = 15;
    }else if(params->dyrankLevel == 3){
        rank_flag_thres = 5;
    }

    
    int max_sub_iter_inc_factor = 10000;
    int max_sub_iter_ceil = 25000;

    int update_max_sub_iter_counter = 0;

    for (lorads_int k = alm_iter_state->outerIter; k <= params->maxALMIter; k++){
        lorads_int rho_inner_iter = 0;
        double rank_update_alpha = 0.1;
        double rank_update_threshold = 0.005;
        double rank_update_current_ema = 0.0;
        double rank_update_old_ema = 0.0;
        lorads_int rank_update_evaluate_interval = 5;
        lorads_int rank_update_counter = 1;
        lorads_int cur_iter_counter = 1;
        if (update_max_sub_iter_counter >= 2){
            update_max_sub_iter_counter = 0;
            MAX_ALM_SUB_ITER += max_sub_iter_inc_factor;
            MAX_ALM_SUB_ITER = LORADS_MIN(MAX_ALM_SUB_ITER, max_sub_iter_ceil);
        }
        while (difficulty != EASY){
            localIter = 0;
            lorads_int if_break = LUtilUpdateCheckEma(&rank_update_current_ema, &rank_update_old_ema, rho_certificate_val, rank_update_alpha, rank_update_threshold, rank_update_evaluate_interval, &rank_update_counter);
            if (!if_break && !params->highAccMode)
            {
                break;
            }
            if (cur_iter_counter >= MAX_ALM_SUB_ITER)
            {
                update_max_sub_iter_counter += 1;
                break;
            }
            if (rank_flag >= rank_flag_thres && (!is_rank_max) && (k - last_outter_iter_start >= 3))
            {
                break;
            }
            if (rho_certificate_val <= rho_certificate_tol)
            {
                break;
            }
            while (rho_certificate_val - rho_certificate_tol > params->endALMSubTol)
            {
                if (localIter % 300 == 0)
                {
                    clearLBFGS = 0;
                }
                aFunc->LBFGSDirection(params, ASolver, ASolver->lbfgsHis, ASolver->var->gradLp, ASolver->var->uLp, ASolver->var->Grad, ASolver->var->U, clearLBFGS);
                aFunc->LBFGSDirUseGrad(ASolver, ASolver->var->uLp, ASolver->var->gradLp, ASolver->var->U, ASolver->var->Grad);
#ifdef FIX_INI_POINT
                double nrm2U = 0.0;
                for (int iCone = 0; iCone < ASolver->nCones; ++iCone){
                    int tempn = ASolver->var->U[iCone]->nRows * ASolver->var->U[iCone]->rank;
                    double temp = nrm2(&tempn, ASolver->var->U[iCone]->matElem, &incx);
                    nrm2U += temp * temp;
                }
                printf("nrm2U: %.20f\n", sqrt(nrm2U));
#endif
                double *q0 = ASolver->var->M1temp;
                // q0 = b - A(RRt)
                LORADS_MEMCPY(q0, ASolver->rowRHS, double, ASolver->nRows);
                axpy(&(ASolver->nRows), &minusOne, ASolver->var->constrValSum, &incx, q0, &incx);

                double p12[2];
                aFunc->ALMCalq12p12(ASolver, ASolver->var->rLp, ASolver->var->uLp, ASolver->var->R, ASolver->var->U, ASolver->var->ARDSum, ASolver->var->ADDSum, p12);
                lorads_int rootNum = ALMLineSearch(alm_iter_state->rho, ASolver->nRows, ASolver->var->dualVar, p12[0], p12[1], q0, ASolver->var->ARDSum, ASolver->var->ADDSum, &tau);
                if (rootNum == 0){
                    retcode = RET_CODE_NUM_ERR;
                    goto END_ALM;
                }
                if (fabs(tau) < params->endTauTol)
                {
                    printf("update rho:%5.8e since tau is too small.\n", tau);
                    alm_iter_state->innerIter++;
                    localIter++;
                    cur_iter_counter++;
                    clearLBFGS++;
                    goto UpdateRho;
                }
                // y = (-grad)
                // set lbfgs one
                aFunc->setAsNegGrad(ASolver, ASolver->var->gradLp, ASolver->var->Grad);
                // update R
                aFunc->ALMupdateVar(ASolver, ASolver->var->rLp, ASolver->var->uLp, ASolver->var->R, ASolver->var->U, tau);
#ifdef FIX_INI_POINT
      printf("tau:%.20f\n", tau);
#endif
                // update gradient
                lagNormSquare = 0.0;
                // update constrValSum first
                double tauSquare = tau * tau;
                axpy(&(ASolver->nRows), &tau, ASolver->var->ARDSum, &incx, ASolver->var->constrValSum, &incx);
                axpy(&(ASolver->nRows), &(tauSquare), ASolver->var->ADDSum, &incx, ASolver->var->constrValSum, &incx);
                aFunc->ALMCalGrad(ASolver, ASolver->var->rLp, ASolver->var->gradLp, ASolver->var->R, ASolver->var->Grad, &lagNormSquare, alm_iter_state->rho);
                aFunc->setlbfgsHisTwo(ASolver, ASolver->var->gradLp, ASolver->var->uLp, ASolver->var->Grad, ASolver->var->U, tau);

                aFunc->updateDimacsALM(ASolver, ASolver->var->R, ASolver->var->R, ASolver->var->rLp, ASolver->var->rLp);
                alm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
                alm_iter_state->l_inf_primal_infeasibility = alm_iter_state->l_1_primal_infeasibility * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf);

                if ((alm_iter_state->l_inf_primal_infeasibility <= params->phase1Tol) && ((alm_iter_state->primal_dual_gap <= params->phase1Tol) || (!params->highAccMode)))
                {
                    alm_iter_state->outerIter = k;
                    alm_iter_state->innerIter += 1;
                    localIter += 1;
                    cur_iter_counter += 1;
                    clearLBFGS += 1;
                    goto END_ALM;
                }
                double lag_grad_norm = sqrt(lagNormSquare);
                rho_certificate_val =  lag_grad_norm/ (1 + ASolver->cObjNrmInf);
                alm_iter_state->innerIter += 1;
                localIter++;
                cur_iter_counter++;
                clearLBFGS++;
                if (localIter > 800){
                    break;
                }
            }
            LORADSUpdateDualVar(ASolver, alm_iter_state->rho);
            aFunc->ALMCalGrad(ASolver, ASolver->var->rLp, ASolver->var->gradLp, ASolver->var->R, ASolver->var->Grad, &lagNormSquare, alm_iter_state->rho);
            rho_certificate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
            if (localIter <= 20){
                difficulty = EASY;
            }
            else if (20 < localIter && localIter <= 100){
                difficulty = MEDIUM;
                rank_flag += 2;
            }
            else if (100 < localIter && localIter < 400){
                difficulty = HARD;
                rank_flag += 3;
            }
            else if (400 <= localIter){
                difficulty = SUPER;
                rank_flag += 4;
            }
            if (difficulty == EASY){
                rank_flag = 0;
            }
            rho_inner_iter ++;
        }
        UpdateRho:
        do{
            alm_iter_state->rho *= rho_update_factor;
            aFunc->ALMCalGrad(ASolver, ASolver->var->rLp, ASolver->var->gradLp, ASolver->var->R, ASolver->var->Grad, &lagNormSquare, alm_iter_state->rho);
            rho_certificate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
            rho_certificate_tol = rho_certificate / alm_iter_state->rho;
        }while(rho_certificate_tol >= rho_certificate_val);
        if (alm_iter_state->rho >= 5e4 && rho_factor_flag < 4){
            rho_update_factor = sqrt(sqrt(rho_update_factor));
            rho_factor_flag = 4;
        }else if(alm_iter_state->rho >= 5e6 && rho_factor_flag < 6){
            rho_update_factor = sqrt(sqrt(rho_update_factor));
            rho_factor_flag = 6;
        } else if (alm_iter_state->rho >= 5e8 && rho_factor_flag < 8){
            rho_update_factor = sqrt(sqrt(rho_update_factor));
            rho_factor_flag = 8;
        }
        difficulty = HARD;
        clearLBFGS = 0;
        alm_iter_state->outerIter = k;
        if (k % 1 == 0)
        {
            if ((alm_iter_state->l_inf_primal_infeasibility <= params->phase1Tol) && ((alm_iter_state->primal_dual_gap <= params->phase1Tol) || (!params->highAccMode)))
            {
                goto END_ALM;
            }
            aFunc->calObj_alm(ASolver);
            LORADSCalDualObj(ASolver);
            aFunc->updateDimacsALM(ASolver, ASolver->var->R, ASolver->var->R, ASolver->var->rLp, ASolver->var->rLp);
            alm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
            alm_iter_state->primal_objective_value = ASolver->pObjVal;
            alm_iter_state->dual_objective_value = ASolver->dObjVal;
            alm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
            alm_iter_state->l_inf_primal_infeasibility = alm_iter_state->l_1_primal_infeasibility * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf);
            alm_iter_state->l_1_dual_infeasibility = 99; // does not evaluate dual error
            alm_iter_state->l_inf_dual_infeasibility = 99; // does not evaluate dual error
            if (alm_iter_state->primal_dual_gap <= params->phase1Tol * 1e-3 && alm_iter_state->l_1_primal_infeasibility <= params->phase1Tol * 1e-3)
            {
                goto PRINT_AND_EXIT;
            }
            // if (alm_iter_state->l_inf_primal_infeasibility <= params->phase2Tol)
            // {
            //     goto PRINT_AND_EXIT;
            // }
            // if (alm_iter_state->primal_dual_gap <= params->phase2Tol * 1e1 && alm_iter_state->l_1_primal_infeasibility <= params->phase2Tol)
            // {
            //     goto PRINT_AND_EXIT;
            // }
            ALMPrintLog(alm_iter_state, LUtilGetTimeStamp() - ori_start);
            if (LUtilGetTimeStamp() - timeSolveStart >= params->timeSecLimit){
                goto PRINT_AND_EXIT;
            }
        }
        if (rank_flag >= rank_flag_thres && !is_rank_max){
            rank_flag = 0;
            if (k - last_outter_iter_start >= 2 && ASolver->nCones <= 1e+10){
                printf("increase the rank, factor:%f.\n", rank_update_factor);
                is_rank_max = AUG_RANK(ASolver, ASolver->var->rankElem, ASolver->nCones, rank_update_factor);
                alm_iter_state->outerIter = k;
                last_outter_iter_start = alm_iter_state->outerIter;
                goto ALG_START;
            }
        }
    }
    END_ALM:
    aFunc->calObj_alm(ASolver);
    LORADSCalDualObj(ASolver);
    aFunc->updateDimacsALM(ASolver, ASolver->var->R, ASolver->var->R, ASolver->var->rLp, ASolver->var->rLp);
    alm_iter_state->primal_objective_value = ASolver->pObjVal;
    alm_iter_state->dual_objective_value = ASolver->dObjVal;
    alm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
    alm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
    alm_iter_state->l_inf_primal_infeasibility = alm_iter_state->l_1_primal_infeasibility * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf);
    alm_iter_state->l_1_dual_infeasibility = 99; // does not evaluate dual error
    alm_iter_state->l_inf_dual_infeasibility = 99; // does not evaluate dual error
    PRINT_AND_EXIT:
        printf("-----------------------------------------------------------------------\n");
        printf("Exit ALM:\n");
        ALMPrintLog(alm_iter_state, LUtilGetTimeStamp() - ori_start);
        printf("-----------------------------------------------------------------------\n");
        return retcode;
}



extern void LORADSCalObjRR_ALM(lorads_solver *ASolver){
    ASolver->pObjVal = 0.0;
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        lorads_sdp_dense *R = ASolver->var->R[iCone];
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        LORADSUVt(ACone->sdp_obj_sum, R, R);
        ACone->objAUV(ACone->coneData, R, R, &ASolver->pObjVal, ACone->sdp_obj_sum);
    }
    ASolver->pObjVal /= ASolver->scaleObjHis;
}

extern void LORADSCalObjRR_ALM_LP(lorads_solver *ASolver){
    ASolver->pObjVal = 0.0;
    ASolver->lpCone->objAUV(ASolver->lpCone->coneData, ASolver->var->rLp, ASolver->var->rLp, &ASolver->pObjVal);
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        lorads_sdp_dense *R = ASolver->var->R[iCone];
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        // can optimize here, zero, dense, sparse difference, really need to UVt for three types matrices?
        LORADSUVt(ACone->sdp_obj_sum, R, R);
        ACone->objAUV(ACone->coneData, R, R, &ASolver->pObjVal, ACone->sdp_obj_sum);
    }
    ASolver->pObjVal /= ASolver->scaleObjHis;
}