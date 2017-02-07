# cs-with-prior-information
This repository provides a Matlab implementation of the ADMM-based solvers for L1-L1 and L1-L2 minimization.
This software provides a Matlab implementation of solvers for the L1-L1 and L1-L2 
minimization problems, as described in

[1] J. F. C. Mota, N. Deligiannis, M. R. D. Rodrigues
    "Compressed Sensing with Prior Information: Optimal Strategies, Geometry, and
    Bounds"
    submitted to IEEE Transactions on Information Theory, 2014
    arXiv: http://arxiv.org/abs/1408.5250

and 

[2] J. F. C. Mota, N. Deligiannis, M. R. D. Rodrigues
    "Compressed Sensing with Side Information: Geometrical Interpretation and 
    Performance Bounds"
    IEEE Global Conf. on Signal and Information Processing (GlobalSIP), pp. 
    512-516, 2014

It also provides code to reproduce the Figures in [1], and a solver for 
Modified-CS, a problem proposed in

[3] N. Vaswani, W. Lu
    "Modified-CS: Modifying Compressive Sensing for Problems With Partially Known 
    Support"
    IEEE Transactions on Signal Processing, Vol. 58, No. 9, 2010

------------
Organization
------------

solvers: contains code for 
           * L1-L1 minimization
           * L1-L2 minimization
           * Modified-CS
         Each folder contains a detailed derivation of the implemented algorithm

createFigures: contains code to reproduce the figures in [1].


If you use this code, please cite [1] and/or [2].


-------
License
-------

Copyright (C) 2017  Joao Mota

This program is free software: you can redistribute it and/or modify it under the 
terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this
program.  If not, see <http://www.gnu.org/licenses/>.

