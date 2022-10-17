# PAPVFCC-HD# PAPVFCC-HM
A Projection Algorithm for Pseudomonotone Vector Fields with Convex Constraints on Hadamard Manifolds

This package includes Matlab codes for a geometric projection algorithm.
They are efficient algorithms for approximating a zero of a pseudomonotone tangent vector field with convex constraints Hadamard manifold. 
The test cases and scripts for running the experiments in paper "A Projection Algorithm for Pseudomonotone Vector Fields with Convex Constraints on
Hadamard Manifolds" by Zhi Zhao, Qin Zeng, Yu-Nong Xu, Ya-Guan Qian, and Teng-Teng Yao, are also included.


1. Main algorithms

HHTWO1.m -- comparition of Algorithm 2.1 with the extragrdient method for Example 4.1

HHTWO2.m --  comparition of Algorithm 2.1 with Korpelevich's method for Example 4.2

2. Test data and codes

EXP.m -- the exponential mapping on H^2

EXPINV-- the inverse of the exponential mapping on H^2 

PHI.m-- backtraking line search on H^2 for Algorithm 2.1 and the the extragrdient method

PHI2.m-- backtraking line search for Korpelevich's method

LORENTZ.m -- the Riemannian metric on H^2 

FF.m -- the value of the underlying vector field of Example 4.1 and Example 4.2 
 
DIST-- the Riemannian distance function on H^2 

PC- the metric projection onto C
