#include "task.hh"

#include <iostream>

#include <typed-geometry/tg.hh>

#include <Eigen/Sparse>
#include <Eigen/SparseCore>

namespace task
{
void init_texture_coordinates(pm::vertex_attribute<tg::pos3> const& position, pm::vertex_attribute<tg::pos2>& texture_coordintate)
{
    auto const& m = position.mesh();

    std::vector<pm::vertex_handle> loop;
    pm::halfedge_handle start_halfedge;
    for (auto const h : m.halfedges())
        if (h.is_boundary())
        {
            start_halfedge = h;
            break;
        }
    TG_ASSERT(start_halfedge.is_valid() && "no boundary found!");

    auto current = start_halfedge;
    do
    {
        loop.push_back(current.vertex_to());
        current = current.next();
    } while (current != start_halfedge);

    // INSERT CODE:
    // 1. Map interior vertices to circle's center (0.5,0.5)
    //
    //    Use the texture_coordintate[pm::vertex_handle] = tg::pos2(<UCoord>,<VCoord>);
    //
    // 2. Map boundary vertices onto a circle with radius 0.5 in texture space (preserve edge length ratio)
    //
    // Hints:
    // - The vector loop contains all vertex handles of the boundary loop in correct order.
    //================================================================================
    //--- start strip ---
    //================================================================================

    for (auto vertex : m.vertices()){
        if (vertex.is_boundary()){ continue; }
        texture_coordintate[vertex] = tg::pos2(0.5, 0.5);
    }

    /*float step_size = (2.f*3.14159)/loop.size();

    for (int i = 0; i < loop.size(); ++i){
        texture_coordintate[loop[i]] = tg::pos2(0.5+0.5*cos(step_size*i), 0.5+0.5*sin(step_size*i));
    }*/

    float total_length = 0;
    std::vector<float> edge_lengths;

    for (int i = 0; i < loop.size()-1; ++i){
        total_length += tg::length(position[loop[i]] - position[loop[i+1]]);
        edge_lengths.push_back(total_length);
    }
    total_length += tg::length(position[loop[loop.size()-1]] - position[loop[0]]);
    edge_lengths.push_back(total_length);

    float ratio = (2*3.14159)/total_length;

    for (int i = 0; i < loop.size(); ++i){
        texture_coordintate[loop[i]] = tg::pos2(0.5+0.5*cos(ratio*edge_lengths[i]), 0.5+0.5*sin(ratio*edge_lengths[i]));
    }

    //================================================================
    //--- end strip ---
    //================================================================
}

void smooth_texcoords(pm::Mesh const& m, int iterations, pm::edge_attribute<float> const& weight, pm::vertex_attribute<tg::pos2>& texture_coordinate)
{
    auto new_texture_coordinate = texture_coordinate;
    for (auto i = 0; i < iterations; ++i)
    {
        for (auto v : m.vertices())
        {
            // INSERT CODE:
            // Iteratively solve equation system by computing the laplace vector of
            // the texture coordinates and adding it onto the current texture coordinates.
            //
            // Hints:
            // - Edge weights can be retrieved via weight[pm::edge_handle]
            // - Parametrization coords for vertices can be retrieved
            //   via texture_coordinate[pm::vertex_handle] as tg::pos2
            // - The new Parametrization coords should be written to
            //   new_texture_coordinate[pm::vertex_handle] as tg::pos2
            //================================================================================
            //--- start strip ---
            //================================================================================

            if (v.is_boundary()){
                continue;
            }

            auto laplacian = tg::pos2::zero;
            float sum_weights = 0;

            for (auto halfedge : v.outgoing_halfedges()){
                auto pi = texture_coordinate[halfedge.vertex_from()];
                auto pj = texture_coordinate[halfedge.vertex_to()];
                laplacian += weight[halfedge.edge()]*(pj - pi);
                sum_weights += weight[halfedge.edge()];
            }
            laplacian /= sum_weights;

            new_texture_coordinate[v] = texture_coordinate[v] + laplacian;

            //================================================================================
            //--- end strip ---
            //================================================================================
        }
        texture_coordinate = new_texture_coordinate;
    }
}

void compute_weights(gp::weight_type type, pm::vertex_attribute<tg::pos3> const& position, pm::edge_attribute<float>& edge_weight)
{
    auto const& m = edge_weight.mesh();

    // Uniform weighting
    if (type == gp::weight_type::uniform)
    {
        edge_weight.clear(1.0);
    }
    // Cotangent weighting
    else if (type == gp::weight_type::cotangent)
    {
        for (auto e : m.edges())
        {
            float w = 0.0f;

            auto h0 = e.halfedgeA();
            auto v0 = h0.vertex_to();
            auto p0 = position(v0);

            auto h1 = e.halfedgeB();
            auto v1 = h1.vertex_to();
            auto p1 = position(v1);

            auto h2 = h0.next();
            auto v2 = h2.vertex_to();
            auto p2 = position(v2);

            auto d0 = normalize(p0 - p2);
            auto d1 = normalize(p1 - p2);
            w += 1.0f / tg::tan(tg::acos(dot(d0, d1)));

            h2 = h1.next();
            p2 = position(h2.vertex_to());
            d0 = normalize(p0 - p2);
            d1 = normalize(p1 - p2);
            w += 1.0f / tg::tan(tg::acos(dot(d0, d1)));

            edge_weight[e] = w;
        }
    }
}

void add_row_to_system(std::vector<Eigen::Triplet<float>>& triplets,
                       pm::vertex_attribute<int> const& sysid,
                       pm::edge_attribute<float> const& weight,
                       pm::vertex_attribute<tg::pos2> const& texture_coordinate,
                       Eigen::VectorXf& rhsu,
                       Eigen::VectorXf& rhsv,
                       pm::vertex_handle origvh)
{
    // INSERT CODE:
    // todo: setup one row of the equation system by pushing back (triplets.push_back(...)) the triplets
    // (of non-zero entries) of the Laplacian for vertex origvh
    //
    // For constrained (boundary) neighbors also add the corresponding right hand side entries to rhsu and rhsv
    //
    // - Use sysid[vh] to get the corresponding row and columns indices in the system matrix
    // - Use texture_coordinate[vh] to retrieve Texture coordinates for rhs
    // - Use weight[eh] to retrive the edge weight
    // - rhsu and rhsv are n-dimensional vectors. You can assign to the i-th entry using rhsv[i] = ...
    //
    // Example for setting the diagonal to a value:
    // triplets.push_back(Eigen::Triplet<float>( sysid[origvh], sysid[origvh], <value>));
    //================================================================================
    //--- start strip ---
    //================================================================================
    float sum_weights = 0;

    for (auto halfedge : origvh.outgoing_halfedges()){
        auto pj = halfedge.vertex_to();
        sum_weights += weight[halfedge.edge()];
        if (pj.is_boundary()){
            rhsu[sysid[origvh]] -= weight[halfedge.edge()] * texture_coordinate[pj].x;
            rhsv[sysid[origvh]] -= weight[halfedge.edge()] * texture_coordinate[pj].y;
            continue;
        }
        triplets.push_back(Eigen::Triplet<float>(sysid[origvh], sysid[pj], weight[halfedge.edge()]));
    }

    triplets.push_back(Eigen::Triplet<float>(sysid[origvh], sysid[origvh], -sum_weights));

    //================================================================================
    //--- end strip ---
    //================================================================================
}

void direct_solve(pm::vertex_attribute<tg::pos3> const& position, pm::edge_attribute<float>& edge_weight, gp::weight_type type, pm::vertex_attribute<tg::pos2>& texture_coordinate)
{
    auto const& m = position.mesh();

    // make sure the boundary has been mapped to a circle (we need these texcoords in add_row_to_system_matrix)
    init_texture_coordinates(position, texture_coordinate);

    // also make sure the weights have been computed
    compute_weights(type, position, edge_weight);

    auto sysid = m.vertices().make_attribute<int>();
    int n_boundary = 0;
    int n_inner = 0;
    for (auto const v : m.vertices())
    {
        if (v.is_boundary())
            ++n_boundary;
        else
            sysid[v] = n_inner++;
    }


    // system matrix
    Eigen::SparseMatrix<float> A(n_inner, n_inner);

    // right hand sides for u and v coordinates
    Eigen::VectorXf rhsu(n_inner);
    rhsu.setZero();
    Eigen::VectorXf rhsv(n_inner);
    rhsv.setZero();

    // resulting texture coordinates for u and v
    Eigen::VectorXf resu(n_inner);
    resu.setZero();
    Eigen::VectorXf resv(n_inner);
    resv.setZero();

    // the matrix is built/initialized from a set of triplets (i.e. (rowid, colid, value))
    /* Info: add_row_to_system should add entries to the vector of triplets, the matrix
     * is build from these triplets below.
     * (see http://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title3 for
     * more details on triplets and setting up sparse matrices)
     */
    std::vector<Eigen::Triplet<float>> triplets;


    // INSERT CODE:
    // for all inner vertices, setup the corresponding row of the linear systems (u and v)
    // using the function add_row_to_system for all inner vertices.
    //================================================================================
    //--- start strip ---
    //================================================================================

    for (auto vertex : m.vertices()){
        if (vertex.is_boundary()){ continue; }
        add_row_to_system(triplets, sysid, edge_weight, texture_coordinate, rhsu, rhsv, vertex);
    }

    //================================================================================
    //--- end strip ---
    //================================================================================

    std::cout << " number of triplets (i.e. number of non-zeros) " << triplets.size() << ", per row " << (triplets.size() / n_inner) << std::endl;

    // now we have all triplets to setup the matrix A
    A.setFromTriplets(triplets.begin(), triplets.end());

    // now we can solve for u and v
    Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> bicg(A); // performs a Biconjugate gradient stabilized method
    resu = bicg.solve(rhsu);
    if (bicg.info() != Eigen::Success)
        std::cerr << "solve failed!" << std::endl;
    resv = bicg.solve(rhsv);
    if (bicg.info() != Eigen::Success)
        std::cerr << "solve failed!" << std::endl;

    // write back to texcoord
    for (auto v : m.vertices())
        if (!v.is_boundary()) // constrained -> skip
            texture_coordinate[v] = tg::pos2(resu[sysid[v]], resv[sysid[v]]);
}

}
