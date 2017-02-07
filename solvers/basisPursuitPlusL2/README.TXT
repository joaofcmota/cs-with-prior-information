
This folder contains two solvers for the problem

   minimize   ||x||_1 + 0.5*beta*||x-w||^2
      x
   subject to A*x = b

where ||.||_1 is the L1-norm and ||.|| is the L2-norm. w is known and is considered side information.


The solver in basisPursuitPlusL2.m is similar to basisPursuitPlusL1 and uses ADMM as the underlying
algorithm. The solver in basisPursuitPlusL2_BB.m applies the Barzilai-Borwein algorithm described in

   E. Birgin, J. Martinez, and M. Raydan, "Nonmonotone Spectral Projected Gradient Methods on Convex
   Sets," SIAM J. Optim., Vol. 10, No. 4, pp. 1196-1211, 2000

According to our tests, basisPursuitPlusL2_BB is faster and more accurate than basisPursuitPlusL2.
Besides, it seems to give the correct solution all the times. In contrast, basisPursuitPlusL2
sometimes fails to converge to the correct solution.

***********************************************
WE RECOMMEND USING basisPursuitPlusL2_BB.m
***********************************************



