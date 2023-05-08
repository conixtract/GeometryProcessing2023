#include "task.hh"
#include <iostream>

#include <typed-geometry/tg.hh>

#define LOG(x) std::cout << x << std::endl;

tg::dir3 task::compute_normal(std::vector<pm::vertex_handle> const& vs, pm::vertex_attribute<tg::pos3> const& position)
{
    // the normal to be computed
    tg::dir3 normal;

    /*
     * INSERT YOUR OWN CODE BETWEEN THE HORIZONTAL LINES BELOW.
     * DO NOT CHANGE CODE ANYWHERE ELSE.
     *
     * Note that there is another function below this function
     * where you have to insert code as well.
     *
     * This function should compute a regression plane for the
     * supplied sequence of points and return an arbitrary one
     * of the two normal vectors of that plane. (The orientation
     * doesn't matter at this point.)
     *
     * Hints:
     *   use tg::mat3 as 3x3 (column major) matrix representation
     *   tg::eigen_decomposition_symmetric may save you a lot of time
     *
     */
    // ----- %< -------------------------------------------------------


    //find center of gravity
    tg::pos3 cog(0, 0, 0);
    for(auto vertex : vs){
        cog += position[vertex];
    }
    cog = cog/vs.size();

    auto positions_on_cog = position;

    LOG(vs.size())
    LOG(position.size())

    //move vertices by cog
    for(auto vertex : vs){
        positions_on_cog[vertex] = tg::pos3(positions_on_cog[vertex] - cog);
    }

    float sum_xx = 0;
    float sum_xy = 0;
    float sum_xz = 0;
    float sum_yy = 0;
    float sum_yz = 0;
    float sum_zz = 0;

    //calculate sums for inertia tensor
    for(auto vertex : vs){
        sum_xx += positions_on_cog(vertex).x * positions_on_cog(vertex).x;
        sum_xy += positions_on_cog(vertex).x * positions_on_cog(vertex).y;
        sum_xz += positions_on_cog(vertex).x * positions_on_cog(vertex).z;
        sum_yy += positions_on_cog(vertex).y * positions_on_cog(vertex).y;
        sum_yz += positions_on_cog(vertex).y * positions_on_cog(vertex).z;
        sum_zz += positions_on_cog(vertex).z * positions_on_cog(vertex).z;
    }

    //compute inertia tensor
    tg::mat3 inertia_tensor;
    inertia_tensor[0] = {sum_xx, sum_xy, sum_xz};
    inertia_tensor[1] = {sum_xy, sum_yy, sum_yz};
    inertia_tensor[2] = {sum_xz, sum_yz, sum_zz};

    auto eigenbasis = tg::eigen_decomposition_symmetric(inertia_tensor);

    //find eigenvector of smallest eigenvalue set result
    if (eigenbasis[0].eigenvalue < eigenbasis[1].eigenvalue){
        normal = tg::dir3(eigenbasis[0].eigenvector);
    }else{ normal = tg::dir3(eigenbasis[1].eigenvector);}

    // ----- %< -------------------------------------------------------
    /*
     *
     * Insert your own code above.
     * NO CHANGES BEYOND THIS POINT!
     *
     */

    return normal;
}

float task::compute_mst_weight(pm::vertex_handle v0, pm::vertex_handle v1, pm::vertex_attribute<tg::pos3> const& position, pm::vertex_attribute<tg::dir3> const& normal)
{
    // this is the weight that you should overwrite
    float weight;

    /*
     * INSERT YOUR OWN CODE BETWEEN THE HORIZONTAL LINES BELOW.
     * DO NOT CHANGE CODE ANYWHERE ELSE.
     *
     * This function should compute the weight of the edge between
     * the two supplied vertices for the purpose of the normal
     * orientation propagation algorithm using the minimum spanning tree.
     *
     */
    // ----- %< -------------------------------------------------------


    // ----- %< -------------------------------------------------------
    /*
     *
     * Insert your own code above.
     * NO CHANGES BEYOND THIS POINT!
     *
     */

    return weight;
}
