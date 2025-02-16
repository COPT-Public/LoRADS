#ifndef LORADS_ELEMENTS_H
#define LORADS_ELEMENTS_H






extern void addDense(double *alpha, void *constrVal, double *vec);

extern void addSparse(double *alpha, void *constrVal, double *vec);

extern void zeroDense(void *constrVal);

extern void zeroSparse(void *constrVal);







#endif