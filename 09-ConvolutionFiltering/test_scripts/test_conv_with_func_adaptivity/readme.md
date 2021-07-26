This folder provides the test scripts for the idea

**Functional Adaptivity** + **L2 projection** + **Convolution postprocesing**

Steps:
1. Build a corase mesh TH
2. Do functinal adaptivity (n times), obtain uh and vh on Th
3. L2 project uh and vh to Th, obtain uH and vH
4. Convolution postprocesing to get uH* and vH*
5. Evaluate the error of uH* and vH*
6. Compare with uH* and vH* obtained simply from TH
