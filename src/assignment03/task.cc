#include "task.hh"

#include <typed-geometry/tg.hh>
#include <iostream>
#define LOG(x) std::cout << x << std::endl;

namespace task
{
pm::edge_attribute<float> compute_weights(pm::Mesh& mesh, pm::vertex_attribute<tg::pos3>& position, bool cotan_weights)
{
    auto weights = mesh.edges().make_attribute<float>();

    // Uniform weighting
    for (auto eh : mesh.edges())
        weights(eh) = 1.0f;

    if (cotan_weights) // Cotangent weighting
    {
        for (auto eh : mesh.edges())
        {
            if (eh.is_boundary())
                continue;

            // INSERT CODE:
            // Compute the cotan weights and store them in the weights attribute
            //--- start strip ---
            // eh is the edge A,B
            //     C
            //    / \
            //   /   \
            //  A-----B
            //   \   /
            //    \ /
            //     D

            auto A = eh.vertexA();
            auto B = eh.vertexB();
            auto C = eh.halfedgeB().next().vertex_to();
            auto D = eh.halfedgeA().next().vertex_to();

            auto vec_C_A = tg::normalize(position(A)-position(C));
            auto vec_C_B = tg::normalize(position(B)-position(C));
            auto vec_D_A = tg::normalize(position(A)-position(D));
            auto vec_D_B = tg::normalize(position(B)-position(D));

            float alpha = acos(tg::dot(vec_D_A, vec_D_B));
            float beta = acos(tg::dot(vec_D_B, vec_D_A));

            weights(eh) = 0.5f * (1.f/tan(alpha) + 1.f/tan(beta));

            //--- end strip ---
        }
    }

    return weights;
}


void compute_new_positions(pm::Mesh& mesh,
                           pm::vertex_attribute<tg::pos3>& position,
                           pm::edge_attribute<float> const& edge_weight,
                           const pm::vertex_attribute<bool>& locked,
                           bool simple_laplace,
                           int iterations)
{
    // Compute new positions using Laplace or Laplace^2 smoothing

    auto new_position = mesh.vertices().make_attribute<tg::pos3>();

    for (int i = 0; i < iterations; ++i)
    {
        // Laplace
        if (simple_laplace)
        {
            for (auto vh : mesh.vertices())
            {
                auto u = tg::vec3::zero;

                // INSERT CODE:
                // Compute the Laplace vector and store the updated position in new_position:
                // new_position(v) = position(v) + 0.5* Laplace(v)
                //--- start strip ---

                auto laplace = tg::vec3::zero;
                for (auto halfedge : vh.outgoing_halfedges()){
                    auto const w = edge_weight(halfedge.edge());
                    auto const pj = position(halfedge.vertex_to());
                    auto const pi = position(vh);

                    laplace += w*(pj - pi);
                }
                laplace /= vh.outgoing_halfedges().size();

                new_position(vh) = position(vh) + 0.5f * laplace;

                //--- end strip ---
            }
        }
        else // bilaplacian smoothing
        {
            // INSERT CODE:
            // Compute the squared Laplacian update
            // 1st: compute Laplaces of positions
            // 2nd: compute Laplaces of Laplacian vectors of all one-ring neighbors
            // 3rd: store updated positions in new_position (use damping factor 0.25 for stability)
            //--- start strip ---

            //--- end strip ---
        }

        // set new positions
        for (auto vh : mesh.vertices())
            if (!locked(vh))
                position(vh) = new_position(vh);
    }
}

}
