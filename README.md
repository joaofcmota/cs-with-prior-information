# Compressed Sensing with Prior Information: Optimal Strategies, Geometry, and Bounds

Matlab solvers for L1-L1 minimization:

<div style="text-align: center;">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{array}[t]{ll}&space;\underset{x}{\text{minimize}}&space;&&space;\|x\|_1&space;&plus;&space;\|x&space;-&space;\overline{x}\|_1&space;\\&space;\text{subject&space;to}&space;&&space;Ax&space;=&space;b,&space;\end{array}" target="_blank"> <img src="https://latex.codecogs.com/gif.latex?\begin{array}[t]{ll}&space;\underset{x}{\text{minimize}}&space;&&space;\|x\|_1&space;&plus;&space;\beta&space;\|x&space;-&space;\overline{x}\|_1&space;\\&space;\text{subject&space;to}&space;&&space;Ax&space;=&space;b&space;\end{array}" title="\begin{array}[t]{ll} \underset{x}{\text{minimize}} & \|x\|_1 + \beta\|x - \overline{x}\|_1 \\ \text{subject to} & Ax = b, \end{array}" /></a>
</div>

L1-L2 minimization:

<div style="text-align: center;">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{array}[t]{ll}&space;\underset{x}{\text{minimize}}&space;&&space;\|x\|_1&space;&plus;&space;\|x&space;-&space;\overline{x}\|_2^2&space;\\&space;\text{subject&space;to}&space;&&space;Ax&space;=&space;b,&space;\end{array}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{array}[t]{ll}&space;\underset{x}{\text{minimize}}&space;&&space;\|x\|_1&space;&plus;&space;\beta&space;\|x&space;-&space;\overline{x}\|_1&space;\\&space;\text{subject&space;to}&space;&&space;Ax&space;=&space;b&space;\end{array}" title="\begin{array}[t]{ll} \underset{x}{\text{minimize}} & \|x\|_1 + \beta\|x - \overline{x}\|_2^2 \\ \text{subject to} & Ax = b, \end{array}" /></a>
</div>

and also for Modified-CS:

<div style="text-align: center;">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{array}[t]{ll}&space;\underset{x}{\text{minimize}}&space;&&space;\|x_{T^c}\|_1\\&space;\text{subject&space;to}&space;&&space;Ax&space;=&space;b\,.&space;\end{array}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{array}[t]{ll}&space;\underset{x}{\text{minimize}}&space;&&space;\|x_{T^c}\|_1\\&space;\text{subject&space;to}&space;&&space;Ax&space;=&space;b\,,&space;\end{array}" title="\begin{array}[t]{ll} \underset{x}{\text{minimize}} & \|x_{T^c}\|_1\\ \text{subject to} & Ax = b\,. \end{array}" /></a>
</div>

The first two problems are analyzed in

1. **[Compressed Sensing with Prior Information: Strategies, Geometry, and Bounds](
    https://doi.org/10.1109/TIT.2017.2695614)**.  
    J. F. C. Mota, N. Deligiannis, M. R. D. Rodrigues.  
    IEEE Transactions on Information Theory, Vol. 63, No. 7, pp. 4472-4496, 2017.  
    [link](https://doi.org/10.1109/TIT.2017.2695614), 
    [arXiv](http://arxiv.org/abs/1408.5250)

and 

2. **[Compressed Sensing with Side Information: Geometrical Interpretation and Performance Bounds](
    http://dx.doi.org/10.1109/GlobalSIP.2014.7032170 )**.  
  **J. F. C. Mota**, N. Deligiannis, M. Rodrigues.  
  IEEE Global Conference on Signal and Information Processing (GlobalSIP),
  Atlanta, 2014. 
  *Session: Information Processing for Big Data.*   
  [link]( http://dx.doi.org/10.1109/GlobalSIP.2014.7032170 ), 
  [arXiv]( http://arxiv.org/abs/1410.2724 )

and Modified-CS was proposed in

3. **[Modified-CS: Modifying Compressive Sensing for Problems With Partially Known Support](
    https://ieeexplore.ieee.org/abstract/document/5471173/)**.  
    N. Vaswani, W. Lu.  
    IEEE Transactions on Signal Processing, Vol. 58, No. 9, 2010.  
    [link](https://ieeexplore.ieee.org/abstract/document/5471173/)

This page also contains code to reproduce the experiments and figures in [1].

## Organization

* solvers: 
  code for L1-L1 minimization,  L1-L2 minimization, Modified-CS. 
  Each folder contains a detailed derivation of the respective algorithm.

* createFigures: code for reproducing the figures in [1].

---

License: [ GPLv3 ]( https://www.gnu.org/licenses/gpl-3.0.en.html )

