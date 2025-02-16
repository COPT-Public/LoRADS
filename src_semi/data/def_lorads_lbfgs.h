#ifndef DEF_LORADS_LBFGS
#define DEF_LORADS_LBFGS


struct lbfgs_node_internal
{
    lorads_int allElem;

    double *s;    // difference of R, Rk - Rk-1
    double *y;    // difference of grad
    double beta;  // 1 / <y, s>
    double alpha; // beta <s, q>
    struct lbfgs_node_internal *next;
    struct lbfgs_node_internal *prev;
};

typedef struct lbfgs_node_internal lbfgs_node;


#endif