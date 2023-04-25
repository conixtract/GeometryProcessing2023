#include "task.hh"

#include <iostream>
#include <queue>

#include <typed-geometry/feature/matrix.hh>
#include <typed-geometry/feature/std-interop.hh>
#include <typed-geometry/feature/vector.hh>

bool task::is_delaunay(polymesh::edge_handle edge, pm::vertex_attribute<tg::pos2> const& position)
{
    auto const va = edge.halfedgeA().next().vertex_to();
    auto const vc = edge.halfedgeA().vertex_to();
    auto const vd = edge.halfedgeB().next().vertex_to();
    auto const vb = edge.halfedgeB().vertex_to();

    /* These are the four points of the triangles incident to edge _eh
            a
           / \
          /   \
        b ----- c
          \   /
           \ /
            d
    */
    tg::pos2 const& a = position[va];
    tg::pos2 const& b = position[vb];
    tg::pos2 const& c = position[vc];
    tg::pos2 const& d = position[vd];

    bool result = true;

    // IMPORTANT: DO NOT ADD ANY CODE OUTSIDE OF THE MARKED CODE STRIPS
    // INSERT CODE:
    // is the edge delaunay or not?
    // -> circum-circle test of the four points (a,b,c,d) OR check if the projected paraboloid is convex
    //--- start strip ---


    //--- end strip ---

    return result;
}

polymesh::vertex_index task::insert_vertex(polymesh::Mesh& mesh, pm::vertex_attribute<tg::pos2>& position, tg::pos2 const& vertex_position, polymesh::face_handle face)
{
    // add vertex and assign it its position
    auto const v = mesh.vertices().add();
    position[v] = vertex_position;

    if (face.is_valid())
    {
        mesh.faces().split(face, v);
        std::cout << "[delaunay] 1:3 Split: vertex " << v.idx.value << " at position " << vertex_position << " inside triangle " << face.idx.value << std::endl;
    }
    else
        return {};

    // IMPORTANT: DO NOT ADD ANY CODE OUTSIDE OF THE MARKED CODE STRIPS
    // INSERT CODE:
    // re-establish Delaunay property
    // ... find edges opposite to the inserted vertex
    // ... are these edges ok? otherwise: flip'em (use mesh.edges().flip(the_edge_to_flip))
    // ... propagate if necessary
    // Hint:
    //   Use is_delaunay(...) (see above) to check the delaunay criterion.
    //   Do not check boundary edges as they do not neighbor two triangles.
    //   You can check if an edge e is boundary by calling e.is_boundary()
    //   You can use an std::queue as a container for edges
    //--- start strip ---


    //--- end strip ---

    return v;
}
