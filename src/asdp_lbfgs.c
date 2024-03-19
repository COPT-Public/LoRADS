#ifdef HEADERPATH
#include "interface/asdp.h"
#include "interface/asdp_utils.h"
#include "src/def_asdp_linsolver.h"
#include "src/asdp_linsolver.h"
#include "src/dense_opts.h"
#include "src/vec_opts.h"
#include "src/def_asdp_linsolver.h"
#include "src/asdp_debug.h"
#else
#include "asdp.h"
#include "asdp_utils.h"
#include "dense_opts.h"
#include "vec_opts.h"
#include "asdp_direct_linsys.h"
#include "asdp_lbfgs.h"
#include "asdp_debug.h"
#endif

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#include <math.h>
double ASDPnthroot(double base, int n) {
    if (base < 0 && n % 2 == 0) {
        return NAN;
    }
    if (base > 0){
        return pow(base, 1.0 / n);
    }else{
        return -pow(-base, 1.0 / n);
    }
    
}

extern int ASDPcubic_equation(double a, double b, double c, double d, double *res){
    double A = b*b - 3*a*c;
    double B = b*c - 9*a*d;
    double C = c*c - 3*b*d;
    double delta = B*B - 4*A*C;
    res[0] = 0.0;
    res[1] = 0.0;
    res[2] = 0.0;
    if(A == 0 && B == 0){
        // one root
        res[0] = ASDP_MAX(res[0], -c/b);
        return 1;
    }else if (delta > 0){
        double Y1 = A*b + 1.5*a*(-B + sqrt(delta));
        double Y2 = A*b + 1.5*a*(-B - sqrt(delta));
        // double Y1_3 = Y1/(2*A*sqrt(A));
        double Y1_3 = ASDPnthroot(Y1, 3);
        double Y2_3 = ASDPnthroot(Y2, 3);
        res[0] = ASDP_MAX(res[0], (-b - Y1_3 - Y2_3)/3/a);
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
