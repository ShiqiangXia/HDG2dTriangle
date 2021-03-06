# HDG2DTriangular 

This project is to implement HDG method to solve PDE problems as well as approximate
output functionals on 2D triangle meshes.

## Features:

1. general 2D domain and triangle meshes

2. Arbitrary precision in Matlab (need the [Multi-precison toolbox](https://www.advanpix.com/))

3. HDG method with Static condensation (source problem and eigenvalue problem)

4. Output functional approximation with the adjoint-based method

5. Adaptive mesh refinement

6. Combination of adaptive strategy with congnvolution filtering 

## To test the code

### 1. Adaptive  HDG

Use **test_functional_adaptive_refinement.m** in 00-TestScripts.

    * Triangle mesh ---> HDG 
    * Error Approximation using HDG local postprocessing
    * Refine mesh --> repeat

### 2. Adaptive HDG + Convolution

Use **conv_adap_test_driver.m** in 09-ConvolutionFiltering/test_scripts.

    * Triangle mesh --> HDG
    * Project to background square mesh
    * Convolve the inner part 
    * Error approx the outer triangle mesh
    * Apative refine the outer mesh untill error small enough


## Paper References
1. [HDG method](https://www.ams.org/journals/mcom/2010-79-271/S0025-5718-10-02334-3/viewer/)
2. Adjoint-based method and Post-Processing
	1. [Adjoint-Based, Superconvergent Galerkin Approximations of Linear Functionals](https://dl.acm.org/doi/abs/10.1007/s10915-017-0507-7)
	2. [An a priori error analysis of adjoint-based super-convergent Galerkin approximations of linear functionals](https://academic.oup.com/imajna/advance-article-abstract/doi/10.1093/imanum/draa102/6104058)

## Some explaination of the code


00-TestScripts

Example of test scripts to run the code

01-ProblemDriver

    * main: enter point of the code
    * Problem driver scripts 
    * Parameter class object and scripts to set paramters

02-Mesh

    * Mesh class object
    * Build mesh and scripts to get useful mesh info

03-Basis&Quad&Integral

    * Polynomial basis on the reference element
    * Gauss Quadratures
    * Matrix for volume and boundary integral

04-HDG
    
    * HDG method for source problems
    * HDG method for eigenvalue problem
    * Local postprocessing 

05-FUncitonal
    
    * Linear functional approximation
    * eigenvalue approxiamtion

06-ErrorAnalysis
    
    * L2 error for solutions
    * Error for functionals
    * Posterior error estimate for funcitonals

07-Visualizaiton

    * Plot 2D solution on given mesh (not implemented yet)
    * Plot elementwise average value

08-Report
    
    * Report problem info
    * Tables for show results

09-ConvolutionFiltering

    * combine adaptive sstrategy with convolution filtering 
    * the main idea is to project solution from triangle mesh to square mesh
      and convolution smooth part while adaptively refining singular part.  
98-Toolbox
    
    * toolbox used

99-RecordResults

    * log

100-Notes

    * notes used for coding


### Toolbox used

1. [Distmesh](http://persson.berkeley.edu/distmesh/) 
2. [ameshref](https://github.com/aschmidtuulm/ameshref)
3. [cprintf](https://www.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-command-window)
4. [labelpoints](https://www.mathworks.com/matlabcentral/fileexchange/46891-labelpoints?s_tid=srchtitle)
