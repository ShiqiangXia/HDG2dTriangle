# HDG2DTriangular 

This project is to implement HDG method to solve PDE problems as well as approximate
output functionals on 2D triangle meshes.

## Features:

1. general 2D domain and triangle meshes

2. Arbitrary precision in Matlab (need the [Multi-precison toolbox](https://www.advanpix.com/))

3. HDG method with Static condensation (source problem and eigenvalue problem)

4. Output functional approximation with the adjoint-based method

5. Adaptive mesh refinement

## To test the code
Run the test files  in 00-TestScripts. 

## Some explaination
00-TestScripts
Example of test scripts to run the code

01-ProblemDriver
    * main: enter point of the code
    * Problem driver scripts 
    * Parameter class object and scripts to set paramters

02-Mesh
    * Mesh class object
    * 
03-Basis&Quad&Integral
04-HDG
05-FUncitonal
06-ErrorAnalysis
07-Visualizaiton
08-Report
98-Toolbox
99-RecordResults
100-Notes


#### Toolbox used

1. [Distmesh](http://persson.berkeley.edu/distmesh/) 
2. [ameshref](https://github.com/aschmidtuulm/ameshref)
3. [cprintf](https://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-command-window)
4. [labelpoints](https://www.mathworks.com/matlabcentral/fileexchange/46891-labelpoints?s_tid=srchtitle)

### References
1. [HDG method] (https://www.ams.org/journals/mcom/2010-79-271/S0025-5718-10-02334-3/viewer/)
2. Adjoint-based method and Post-Processing
	1. [Adjoint-Based, Superconvergent Galerkin Approximations of Linear Functionals](https://dl.acm.org/doi/abs/10.1007/s10915-017-0507-7)
	2. [An a priori error analysis of adjoint-based super-convergent Galerkin approximations of linear functionals](https://academic.oup.com/imajna/advance-article-abstract/doi/10.1093/imanum/draa102/6104058)