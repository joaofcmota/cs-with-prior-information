
Implements an algorithm to solve the problem proposed in

N. Vaswani, W. Lu, "Modified-CS: Modifying Compressive Sensing for Problems
With Partially Known Support," IEEE Transactions on Signal Processing, Vol.
58, No. 9, September 2010.

That problem is

	minimize    ||x_Tc||_1        (Modified-CS)
       x
	subject to  Ax = b

where Tc denotes the complement of T, a known estimate of the support of x*.

We solve this problem with ADMM. 


