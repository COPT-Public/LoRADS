## **[LoRADS](https://github.com/COPT-Public/LoRADS)**

#### A Low Rank ADMM Splitting Approach for Semidefinite Programming

LoRADS is an enhanced first-order method solver for low rank Semi-definite programming problems (SDPs). LoRADS is written in ANSI C and is maintained by Cardinal Operations COPT development team. More features are still under active development.



#### Optimization Problem:

LoRADS focus on the following problem:

$$
\min_{\mathcal{A} \mathbf{X} = \mathbf{b}, \mathbf{X}\in \mathbb{S}_+^n} \left\langle \mathbf{C}, \mathbf{X} \right\rangle
$$

##### Features of the problem:

- linear objective
- affine constraints
- positive-semidefinite variables



#### Current release:

LoRADS is now under active development and release a pre-built binary (v1.0.0) which reads SDPA (.dat-s) format files, solves the SDP problems. Users testing solvers in Linux can download the binary from **[the release site](https://github.com/COPT-Public/LoRADS/releases/tag/v1.0.0)** 



#### Getting started:

After downloading the binary form the release site, you could execute

```sh
unzip LoRADS_v_1_0_0-alpha.zip
chmod +x LoRADS_v_1_0_0-alpha
```

Now,  by running

```sh
LoRADS_v_1_0_0-alpha SDPAFILE.dat-s
```

we can solve SDPs represented in standard SDPA format.



If everything goes well, we would see logs like below:

```
-----------------------------------------------------------
  L         OOO      RRRR       A      DDDD       SSS
  L        O   O     R   R     A A     D   D     S
  L        O   O     RRRR     AAAAA    D   D      SSS
  L        O   O     R  R     A   A    D   D         S
  LLLLL     OOO      R   R    A   A    DDDD       SSS
-----------------------------------------------------------
Input file name: SDPAFILE.dat-s
timesLogRank   : 2.00
phase1Tol      : 1.00e-03
initRho        : 1/n
rhoMax         : 5000.00
rhoFreq        : 5
rhoFactor      : 1.20
heursitcFactor : 1.00
maxIter        : 10000
timeSecLimit   : 10000.00
-----------------------------------------------------------
Reading SDPA file in 0.093176 seconds 
nConstrs = 2964, sdp nBlks = 22, lp Cols = 0
Pre-solver starts 
  Processing the cones 
  End preprocess 
**First using BM method as warm start
Iter:0 objVal:1.13715e+03 dualObj:0.00000e+00 ConstrVio(1):3.01346e+00 ConstrVio(Inf):2.71372e+01 PDGap:9.99121e-01 rho:0.02 minIter:0 trace:6329.42 Time:0.02
Iter:1 objVal:2.56699e+02 dualObj:8.34520e+01 ConstrVio(1):1.15143e+00 ConstrVio(Inf):1.03690e+01 PDGap:5.07830e-01 rho:0.04 minIter:3 trace:5421.49 Time:0.05
Iter:2 objVal:1.26633e+03 dualObj:8.46344e+01 ConstrVio(1):2.60669e-01 ConstrVio(Inf):2.34741e+00 PDGap:8.74058e-01 rho:0.08 minIter:8 trace:5072.61 Time:0.09
...
Iter:8 objVal:6.74243e+01 dualObj:6.73879e+01 ConstrVio(1):1.03085e-03 ConstrVio(Inf):9.28310e-03 PDGap:2.67973e-04 rho:5.18 minIter:275 trace:396.97 Time:1.76
Iter:9 objVal:6.76276e+01 dualObj:6.83311e+01 ConstrVio(1):5.09176e-04 ConstrVio(Inf):4.58529e-03 PDGap:5.13659e-03 rho:10.36 minIter:408 trace:397.59 Time:2.57
**Complete ALM+BM warm start
objVal:6.79500e+01 dualObj:6.81592e+01 ConstrVio:1.10875e-04 Assym:0.00000e+00 DualInfe:1.00000e+00 PDGap:1.52578e-03 rho:10.36 minIter:830 Time:5.12
BM 2 ADMM time consuming :0.000022
**Change method into ADMM Split method
Iter:0 objVal:6.79624e+01 dualObj:6.80730e+01 ConstrVio(1):7.21993e-05 ConstrVio(Inf):6.50177e-04 PDGap:8.07182e-04 rho:12.44 cgIter:1235 trace:398.15 Time:0.21
Iter:1 objVal:6.79678e+01 dualObj:6.80502e+01 ConstrVio(1):3.91495e-05 ConstrVio(Inf):3.52554e-04 PDGap:6.01839e-04 rho:12.44 cgIter:2549 trace:398.16 Time:0.42
Iter:2 objVal:6.79700e+01 dualObj:6.79961e+01 ConstrVio(1):2.43977e-05 ConstrVio(Inf):2.19709e-04 PDGap:1.90828e-04 rho:12.44 cgIter:3877 trace:398.16 Time:0.65
...
Iter:19 objVal:6.79539e+01 dualObj:6.79445e+01 ConstrVio(1):1.60800e-06 ConstrVio(Inf):1.44805e-05 PDGap:6.87056e-05 rho:21.49 cgIter:26245 trace:398.16 Time:4.48
Iter:20 objVal:6.79541e+01 dualObj:6.79417e+01 ConstrVio(1):8.59376e-07 ConstrVio(Inf):7.73895e-06 PDGap:9.01021e-05 rho:25.79 cgIter:27525 trace:398.16 Time:4.69
-----------------------------------------------------------------------
End Program due to reaching `final terminate criteria`:
-----------------------------------------------------------------------
Objective function Value are:
	 1.Primal Objective:            : 6.80e+01
	 2.Dual Objective:              : 6.79e+01
Dimacs Error are:
	 1.Constraint Violation(1)      : 8.59e-07
	 2.Dual Infeasibility(1)        : 1.60e-04
	 3.Primal Dual Gap              : 9.01e-05
	 4.Primal Variable Semidefinite : 0.00e+00
	 5.Constraint Violation(Inf)    : 7.74e-06
	 6.Dual Infeasibility(Inf)      : 3.13e-03
-----------------------------------------------------------------------
Solving SDPAFILE.dat-s in 9.837610 seconds 
Solving + calculate full dual infeasibility in 9.880412 seconds  
```

#### Parameters

LoRADS provides users with customizable parameters to fine-tune the solving process according to specific problem requirements (if needed). Below is a detailed description of each parameter:

| **Parameter**   | **Description**                                                                           | **Type** | **Default Value** |
| --------------- | ----------------------------------------------------------------------------------------- | -------- | ----------------- |
| timesLogRank    | Multiplier for the O(log(m)) rank calculation (rank = **timesLogRank** $\times$ log(m)).  | float    | 2.0               |
| phase1Tol       | Tolerance for ending Phase I.                                                             | float    | 1e-3              |
| initRho         | Initial value for the penalty parameter $\rho$.                                           | float    | 1/n               |
| rhoMax          | Maximum value for the penalty parameter $\rho$.                                           | float    | 5000.0            |
| rhoFreq         | Frequency of increasing $\rho$ (increased every **rhoFreq** iterations).                  | int      | 5                 |
| rhoFactor       | Multiplier for increasing $\rho$ ($\rho =$ **rhoFactor** $\times$ $\rho$).                | float    | 1.2               |
| heuristicFactor | Heuristic factor applied when switching to Phase II.                                      | float    | 1.0               |
| maxIter         | Maximum iteration number for the ADMM algorithm.                                          | int      | 10000             |
| timeSecLimit    | Solving time limitation in seconds.                                                       | float    | 10000.0           |

For example, to set **timesLogRank** to 1.0 and solve a problem, we can execute

```
LoRADS_v_1_0_0-alpha SDPAFILE.dat-s --timesLogRank 1.0
```

#### Contributing

LoRADS is still in its preliminary release and will start accepting pull requests in a future release.


#### Developers

LoRADS is developed by 

- Zhenwei Lin: zhenweilin@163.sufe.edu.cn
- Qiushi Han: joshhan2@illinois.edu


#### Reference

- Qiushi Han, Chenxi Li, Zhenwei Lin, Caihua Chen, Qi Deng, Dongdong Ge, Huikang Liu, and Yinyu Ye. "A Low-Rank ADMM Splitting Approach for Semidefinite Programming." *arXiv preprint arXiv:2403.09133* (2024).
- Burer, Samuel, R. Monteiro, and Changhui Choi. "SDPLR 1.03-beta User’s Guide (short version)(2009)." *URL http://sburer.github.io/files/SDPLR-1.03-beta-usrguide.pdf*.
