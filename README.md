# Numerical-Integration

Numerical Integration in 1D, 2D and 3D using Gauss-Chebychev-Quadrature on a rectangular grid. Function evaluations are done using simplistic C++ lambda functions. 

For a quick theoretical introduction on Gaussian-Quadrature please visit i.e.: https://en.wikipedia.org/wiki/Gaussian_quadrature or see Paul DeVries: Computational Physics Chapter 4

## Usage & Compilation on Linux: 
*Be sure that your system is compatible to C++14 and to have OpenMP installed for parallelization!*
```
make
./MultiDim_Gauss.out
```  
Or compile from source using g++: 
```
g++ MultiDim_Gauss.cpp -fopenmp -std=c++1y
```
  
  
**Any Questions or Bug Reports?** Send me an e-Mail: dominik.lindorfer@jku.at
