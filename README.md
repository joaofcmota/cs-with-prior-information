# Compressed Sensing with Prior Information: Optimal Strategies, Geometry, and Bounds

Matlab solvers for <img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/839a0dc412c4f8670dd1064e0d6d412f.svg?invert_in_darkmode" align=middle width=13.40191379999999pt height=22.831056599999986pt/>-<img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/839a0dc412c4f8670dd1064e0d6d412f.svg?invert_in_darkmode" align=middle width=13.40191379999999pt height=22.831056599999986pt/> minimization:

<p align="center"><img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/e446878229a2485dddbb4998ee40660f.svg?invert_in_darkmode" align=middle width=213.5559888pt height=40.84012515pt/></p>

<img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/839a0dc412c4f8670dd1064e0d6d412f.svg?invert_in_darkmode" align=middle width=13.40191379999999pt height=22.831056599999986pt/>-<img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/336fefe2418749fabf50594e52f7b776.svg?invert_in_darkmode" align=middle width=13.40191379999999pt height=22.831056599999986pt/> minimization:

<p align="center"><img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/5884b1a5a4a2cfacb5e28aa2347864ba.svg?invert_in_darkmode" align=middle width=213.5559888pt height=41.89223445pt/></p>

and also for Modified-CS:

<p align="center"><img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/737e86dfd03ff46ed0f5ad4bb70db5ad.svg?invert_in_darkmode" align=middle width=143.48172465pt height=40.84012515pt/></p>

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
  code for <img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/839a0dc412c4f8670dd1064e0d6d412f.svg?invert_in_darkmode" align=middle width=13.40191379999999pt height=22.831056599999986pt/>-<img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/839a0dc412c4f8670dd1064e0d6d412f.svg?invert_in_darkmode" align=middle width=13.40191379999999pt height=22.831056599999986pt/> minimization, <img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/839a0dc412c4f8670dd1064e0d6d412f.svg?invert_in_darkmode" align=middle width=13.40191379999999pt height=22.831056599999986pt/>-<img src="https://github.com/joaofcmota/cs-with-prior-information/svgs/336fefe2418749fabf50594e52f7b776.svg?invert_in_darkmode" align=middle width=13.40191379999999pt height=22.831056599999986pt/> minimization, Modified-CS. 
  Each folder contains a detailed derivation of the respective algorithm.

* createFigures: code for reproducing the figures in [1].

---

License: [ GPLv3 ]( https://www.gnu.org/licenses/gpl-3.0.en.html )

