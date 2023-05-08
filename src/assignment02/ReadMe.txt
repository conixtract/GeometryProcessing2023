Geometry Processing, Assignment 2
=================================

Given a point cloud, estimate a normal vector for every point. Orientation
of the normals must be consistent.

Hints:
------
tg::covariance_matrix and tg::eigen_decomposition_symmetric can save you
quite a lot of time

In the directory ./data you can find two point clouds: sphere.off and
tetra_thing.off. 

What you should implement:
--------------------------

In task.cc find the *two* functions
compute_normal()
and
compute_mst_weight().
Read *both*, the comments above the respective functions and the comments
inside the functions for instructions on what to implement.

Hand in your solution:
--------------------------
1) Create a zip file containing a MEMBERS.txt and the files you changed 
   in the exercise. For this exercise you should only change and 
   upload task.cc together with the members file.

2) Upload the zip file to the moodle.