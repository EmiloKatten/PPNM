PROJECT 26:
Eigenvalues via rootfinding: Variational method with Lagrange multipliers


Using the analytical Jacobian, rootfinding via Newton's method is greatly facilitated.

In this project:
- I've constructed the optimized rootfinding method. It has been testet on a simple f(x,y) (Himmelblau's) function for debugging.
- The newton method was timed on increasing sizes of symmetric matrices.
- I compared the Newton computation time with that of EVD
- From the EVD homework, I've used my newton rootfinding method to determine the lowest states of hydrogen (s1, s2, s3)
- I perturbed the inital guesses for eigenvalues and eigenvectors when finding the eigenpairs of s2, using s1 and s3.
  This was done to analyse how much Newton rootfinding depends on good initial guesses, and how many calls must be done.
  
The 'Out.txt' file presents and discusses the findings of this project in detail.


