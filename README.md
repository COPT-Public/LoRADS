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

LoRADS is now under active development and release a pre-built binary (v1.0.0) which reads SDPA (.dat-s) format files, solves the SDP problems. Users testing solvers in Linux can download the binary from



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
Input file name: G11.dat-s
timesLogRank   : 2.00
phase1Tol      : 1.00e-03
initRho        : 1/n
rhoMax         : 5000.00
rhoFreq        : 5
rhoFactor      : 1.20
heursitcFactor : 1.00
maxIter        : 100000
timeSecLimit   : 10000.00
-----------------------------------------------------------
Reading SDPA file in 0.003098 seconds 
nConstrs = 800, sdp nBlks = 1, lp Cols = 0

ASDP: software for semi-definite programming 

---------------------------------------------
Pre-solver starts 
  Processing the cones 
  End preprocess 
**Detected MaxCut problem: set phase1Tol -> 1e-2 and heuristicFactor -> 10
**First using BM method as warm start
Iter:0 objVal:-4.79077e+01 dualObj:0.00000e+00 ConstrVio(1):5.51657e-02 ConstrVio(Inf):2.20939e+01 PDGap:9.79553e-01 rho:0.00 minIter:0 trace:1890.10 Time:0.00
Iter:1 objVal:-6.83615e+04 dualObj:-1.39874e+03 ConstrVio(1):1.86971e+00 ConstrVio(Inf):7.48819e+02 PDGap:9.59885e-01 rho:0.01 minIter:24 trace:39351.93 Time:0.01
Iter:2 objVal:-9.11715e+03 dualObj:-1.41716e+03 ConstrVio(1):2.04835e-01 ConstrVio(Inf):8.20362e+01 PDGap:7.30874e-01 rho:0.01 minIter:27 trace:4485.24 Time:0.01
Iter:3 objVal:-4.25715e+02 dualObj:-1.41147e+03 ConstrVio(1):2.63646e-02 ConstrVio(Inf):1.05590e+01 PDGap:5.36265e-01 rho:0.02 minIter:30 trace:230.19 Time:0.01
Iter:4 objVal:-3.68034e+02 dualObj:-1.39932e+03 ConstrVio(1):2.74979e-02 ConstrVio(Inf):1.10129e+01 PDGap:5.83189e-01 rho:0.08 minIter:31 trace:192.64 Time:0.01
Iter:5 objVal:-2.65777e+02 dualObj:-1.34550e+03 ConstrVio(1):3.00399e-02 ConstrVio(Inf):1.20310e+01 PDGap:6.69688e-01 rho:0.16 minIter:33 trace:127.25 Time:0.01
Iter:6 objVal:-7.63862e+02 dualObj:-1.28092e+03 ConstrVio(1):1.96701e-02 ConstrVio(Inf):7.87787e+00 PDGap:2.52742e-01 rho:0.32 minIter:39 trace:396.37 Time:0.01
Iter:7 objVal:-1.31614e+03 dualObj:-1.26977e+03 ConstrVio(1):7.50862e-03 ConstrVio(Inf):3.00720e+00 PDGap:1.79233e-02 rho:0.64 minIter:48 trace:765.17 Time:0.01
Iter:8 objVal:-1.28442e+03 dualObj:-1.25904e+03 ConstrVio(1):2.18708e-03 ConstrVio(Inf):8.75926e-01 PDGap:9.97453e-03 rho:1.28 minIter:78 trace:797.14 Time:0.02
Iter:9 objVal:-1.26149e+03 dualObj:-1.25803e+03 ConstrVio(1):3.72056e-04 ConstrVio(Inf):1.49008e-01 PDGap:1.37379e-03 rho:2.56 minIter:120 trace:799.79 Time:0.03
Iter:10 objVal:-1.25851e+03 dualObj:-1.25819e+03 ConstrVio(1):3.54330e-05 ConstrVio(Inf):1.41909e-02 PDGap:1.27154e-04 rho:5.12 minIter:182 trace:800.03 Time:0.04
**Complete ALM+BM warm start
objVal:-1.25824e+03 dualObj:-1.25819e+03 ConstrVio:7.12325e-06 Assym:0.00000e+00 DualInfe:1.00000e+00 PDGap:1.85166e-05 rho:5.12 minIter:183 Time:0.04
BM 2 ADMM time consuming :0.000002
**Change method into ADMM Split method
Iter:0 objVal:-1.25820e+03 dualObj:-1.25819e+03 ConstrVio(1):1.08384e-06 ConstrVio(Inf):4.34079e-04 PDGap:2.54865e-06 rho:61.44 cgIter:24 trace:800.00 Time:0.00
Iter:1 objVal:-1.25819e+03 dualObj:-1.25819e+03 ConstrVio(1):7.53851e-08 ConstrVio(Inf):3.01917e-05 PDGap:1.52230e-07 rho:61.44 cgIter:53 trace:800.00 Time:0.00
Iter:2 objVal:-1.25819e+03 dualObj:-1.25819e+03 ConstrVio(1):1.57822e-08 ConstrVio(Inf):6.32078e-06 PDGap:2.38178e-08 rho:61.44 cgIter:70 trace:800.00 Time:0.00
-----------------------------------------------------------------------
End Program due to reaching `BM terminate criteria`:
-----------------------------------------------------------------------
Objective function Value are:
	 1.Primal Objective:            : -1.26e+03
	 2.Dual Objective:              : -1.26e+03
Dimacs Error are:
	 1.Constraint Violation(1)      : 1.58e-08
	 2.Dual Infeasibility(1)        : 2.95e-07
	 3.Primal Dual Gap              : 2.38e-08
	 4.Primal Variable Semidefinite : 0.00e+00
	 5.Constraint Violation(Inf)    : 6.32e-06
	 6.Dual Infeasibility(Inf)      : 2.18e-04
-----------------------------------------------------------------------
Solving G11.dat-s in 0.052254 seconds 
Solving + calculate full dual infeasibility in 0.074157 seconds 
```

#### Parameters

LoRADS provides users with customizable parameters to fine-tune the solving process according to specific problem requirements (if needed). Below is a detailed description of each parameter:

| **Parameter**  | **Description**                                                                              | **Type** | **Default Value** |
|----------------|----------------------------------------------------------------------------------------------|----------|-------------------|
| timesLogRank   | Multiplier for the $O(\log(m))$ rank calculation (rank = **timesLogRank** \times $\log(m)$). | float    | 2.0               |
| phase1Tol      | Tolerance for ending Phase I.                                                                | float    | 1e-3              |
| initRho        | Initial value for the penalty parameter $\rho$.                                              | float    | 1/n               |
| rhoMax         | Maximum value for the penalty parameter $\rho$.                                              | float    | 5000.0            |
| rhoFreq        | Frequency of increasing $\rho$ (increased every **rhoFreq** iterations).                     | int      | 5                 |
| rhoFactor      | Multiplier for increasing $\rho$ ( $\rho = $ **rhoFactor** $\times \rho$ ).                  | float    | 1.2               |
| heursitcFactor | Heuristic factor applied when switching to Phase II.                                         | float    | 1.0               |
| maxIter        | Maximum iteration number for the ADMM algorithm.                                             | int      | 10000             |
| timeSecLimit   | Solving time limitation in seconds.                                                          | float    | 10000.0           |

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

- Han, Qiushi, Chenxi Li, Zhenwei Lin, Caihua Chen, Qi Deng, Dongdong Ge, Huikang Liu, and Yinyu Ye. "A Low-Rank ADMM Splitting Approach for Semidefinite Programming." *arXiv preprint arXiv:2403.09133* (2024).
- Burer, Samuel, R. Monteiro, and Changhui Choi. "SDPLR 1.03-beta User’s Guide (short version)(2009)." *URL http://sburer. github. io/files/SDPLR-1.03-beta-usrguide. pdf*.











