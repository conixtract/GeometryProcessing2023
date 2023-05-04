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

    //lamda for calculating the minor of a matrix
    auto minorOfMatrix = [](tg::mat4 matrix, int row, int col){
        tg::mat3 minor;
        int subi = 0, subj = 0;
        for (int i = 0; i < 4; i++) {
            if (i == row-1)
                continue;
            subj = 0;
            for (int j = 0; j < 4; j++) {
                if (j == col-1)
                    continue;
                minor[subi][subj] = matrix[i][j];
                subj++;
            }
            subi++;
        }
        return tg::determinant(minor);
    };

    //center of circle
    tg::pos2 x = {0, 0};

    //setup circle representation matrix
    //source: https://math.stackexchange.com/questions/213658/get-the-equation-of-a-circle-when-given-3-points
    tg::mat<4, 4, float> circle_matrix;

    circle_matrix[0] = {x[0]*x[0] + x[1]*x[1], x[0], x[1], 1};
    circle_matrix[1] = {d[0]*d[0] + d[1]*d[1], d[0], d[1], 1};
    circle_matrix[2] = {b[0]*b[0] + b[1]*b[1], b[0], b[1], 1};
    circle_matrix[3] = {c[0]*c[0] + c[1]*c[1], c[0], c[1], 1};

    x[0] = 0.5f * (minorOfMatrix(circle_matrix, 1, 2)/minorOfMatrix(circle_matrix, 1, 1));
    x[1] = -0.5f * (minorOfMatrix(circle_matrix, 1, 3)/minorOfMatrix(circle_matrix, 1, 1));
    
    auto radius = tg::norm(x - d, 2.f);

    //check empty circumcircle condition
    auto a_distance_to_center = tg::norm(x - a,2.f);

    if (a_distance_to_center < radius){ result = false; }

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
    
    //lambda for checking edges of a given face
    auto checkEdgesOfFace = [](polymesh::face_handle const face, std::queue<polymesh::edge_handle>& edge_queue, pm::vertex_attribute<tg::pos2> const& position){
        for(auto edge_handle : face.edges()){
            if (is_delaunay(edge_handle, position) || edge_handle.is_boundary()){ continue; }
            edge_queue.push(edge_handle);
        }
    };

    std::queue<polymesh::edge_handle> edges_to_flip;

    //find opposite edges and check if delaunay
    for(auto halfedge_handle : v.outgoing_halfedges()){
        auto const opposite_edge = halfedge_handle.next().edge();

        if (is_delaunay(opposite_edge, position) || opposite_edge.is_boundary()){ continue; }
        edges_to_flip.push(opposite_edge);
    }

    //propagate to other edges
    while(!edges_to_flip.empty()){
        auto current_edge = edges_to_flip.front();
        edges_to_flip.pop();
        mesh.edges().flip(current_edge);

        //check edges of face A
        checkEdgesOfFace(current_edge.faceA(), edges_to_flip, position);
        //check edges of face B
        checkEdgesOfFace(current_edge.faceB(), edges_to_flip, position);
    }

    //--- end strip ---

    return v;
}
