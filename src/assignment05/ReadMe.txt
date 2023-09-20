Geometry Processing, Practical Assignment 5 - Parametrization
=============================================================

Parametrization / Texturing


In task.cc you should implement:
-------------------------------------------------------------

1) Map boundary vertices in 3D onto circle in texture space, and
   interior vertices onto circle's center.

2) Smooth the parameterization by Laplacian relaxation,
   compare uniform and cot weighting.

3) Implement the direct (equation system-based) solution by adding
   the missing pieces to direct_solve(...) and add_row_to_system(...).


Hints:

 Let nv be the number of vertices in the mesh and nv_bdry the number of boundary vertices, 
 then the numer of unknowns is (nv-nv_bdry) and the size of the square equation system is (nv-nv_bdry)x(nv-nv_bdry).
 The inner vertices of the mesh need to be re-indexed to the range (0..nv - nv_bdry - 1) when setting up the system. 
 This re-indexing is already setup, use the function sysid[vh] to get the row/column index of the vertex handle vh of an inner vertex.

