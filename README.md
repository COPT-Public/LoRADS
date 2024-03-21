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

LoRADS is now under active development and release a pre-built binary (v1.0.0) which reads SDPA (.dat-s) format files, solves the SDP problem. Users testing solvers in Linux can download the binary from



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



If everything goes well, we would see logs like below.



#### Developers

HDSDP is developed by 

- Zhenwei Lin: zhenweilin@163.sufe.edu.cn
- Qiushi Han: joshhan2@illinois.edu



#### Reference

- Han, Qiushi, Chenxi Li, Zhenwei Lin, Caihua Chen, Qi Deng, Dongdong Ge, Huikang Liu, and Yinyu Ye. "A Low-Rank ADMM Splitting Approach for Semidefinite Programming." *arXiv preprint arXiv:2403.09133* (2024).
- Burer, Samuel, R. Monteiro, and Changhui Choi. "SDPLR 1.03-beta User’s Guide (short version)(2009)." *URL http://sburer. github. io/files/SDPLR-1.03-beta-usrguide. pdf*.











