# HDG2DTriangular 

This project is to implement HDG method on general (unstructured or structured) 2D triangular meshes and with mesh adaptivity.

## To test the code
Run the test files  in 00-TestScripts. 

### Features:

1. general 2D domain and  triangle meshes

2. Arbitrary precision in Matlab (Multi-precison tool box)

3. HDG method with Static condensation (source problem and eigenvalue problem)

4. Adaptive mesh refinement


#### Toolbox used

1. Distmesh 
2. ameshref
3. cprintf
4. labelpoints

### References
1. [HDG method] (https://www.ams.org/journals/mcom/2010-79-271/S0025-5718-10-02334-3/viewer/)
2. Adjoint-based method and Post-Processing
	1. [Adjoint-Based, Superconvergent Galerkin Approximations of Linear Functionals](https://dl.acm.org/doi/abs/10.1007/s10915-017-0507-7)
	2. [An a priori error analysis of adjoint-based super-convergent Galerkin approximations of linear functionals](https://academic.oup.com/imajna/advance-article-abstract/doi/10.1093/imanum/draa102/6104058)