## Welcome to the Riemann Solver Repository

### This code is managed by the Rodriguez Flow Research group at Brown University
### Developer: Mauro Rodriguez Jr. (mauro_rodriguez@brown.edu)

* This code is for educational and testing purposes only
* The code solves the exact Riemann problem with for a prescribed left and right state to a prescribed precision
* The P-star method is evaluated for the exact solution
* The approximate Riemann solver of Roe (Roe solver) with an entropy correction and flux limiter is also computed, 2nd order accuracy is obtained
* Superbee and Minmod flux limiters are available, additional flux limiters can be easily added and tested

#### To do: 
1. Create a more general Riemann solver for a non-constant B (liquid stiffness)
2. Create a multiphase exact Riemann solver
3. Add HLL, HLLC, and HLLE approximate Riemann solvers

