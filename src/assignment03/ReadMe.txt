Geometry Processing, Assignment 3
=================================

Mesh smoothing

What you should implement:
--------------------------
1) Computation of the cotangent weights.

2) Smoothing of the mesh using Laplace operator
   new_position(v) = position(v) + 0.5 * Laplace(v)

3) Smoothing of the mesh using Laplace^2 operator
   new_position(v) = v - 0.25 * Laplace^2(v)
   1st: compute the Laplace of the positions involving a vertex's one-ring
   2nd: compute Laplace of Laplacian vectors of all one-ring neighbors
   3rd: store updated positions in new_position
   damping factor 0.25 ensures convergence
