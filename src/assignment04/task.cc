#include "task.hh"

#include <polymesh/properties.hh>
#include <typed-geometry/tg.hh>
#include "QuadricT.hh"

#include <queue>

namespace task
{
bool is_collapse_legal(pm::vertex_attribute<tg::pos3>& position, pm::face_attribute<tg::vec3> const& normal, pm::halfedge_handle heh, tg::angle32 max_angle)
{
    // collect vertices
    auto const v0 = heh.vertex_from();
    auto const v1 = heh.vertex_to();

    // collect faces that would be removed by the collapse
    auto const fl = heh.face();
    auto const fr = heh.opposite_face();

    // backup point positions
    auto const p0 = position(v0);
    auto const p1 = position(v1);

    // topological test
    if (!pm::can_collapse(heh))
        return false;

    // test boundary stuff
    if (v0.is_boundary() && !v1.is_boundary())
        return false;

    bool collapseOK = true;

    /*
     * INSERT YOUR OWN CODE BETWEEN THE HORIZONTAL LINES BELOW.
     * DO NOT CHANGE CODE ANYWHERE ELSE.
     *
     * Note that there are other functions below this function
     * where you have to insert code as well.
     *
     * This function should simulate a halfedge collapse by setting the
     * position of v0 to the one of v1. The test should check whether
     * the normal vector of the non-degenerate triangles change by
     * more than a given threshold.
     * If this is the case set variable collapseOK to false.
     * Otherwise leave it as true.
     * In any case: undo the simulation in the end. This is important!
     *
     * Note:
     *
     * v0, v1, fl, fr, p0, and p1 are already declared and initialized
     * (see code above)
     *
     * You can:
     *
     * set the position using the "position" attribute
     * access the stored face normal (from before the simulated collapse) using the "normal" attribute.
     * compute a normal for the new position using pm::triangle_normal(fh, position)
     */

    // ----- %< -------------------------------------------------------

    position(v0) = position(v1);

    std::queue<pm::face_handle> faces_to_check;

    for (auto face : v0.faces()){
        if(face != fl && face != fr){ faces_to_check.push(face); }
    }

    for (int i = 0; i < faces_to_check.size(); ++i){
        auto face = faces_to_check.front();
        faces_to_check.pop();

        auto angle_diviation = acos(tg::dot(pm::triangle_normal(face, position), normal(face)));
        if(angle_diviation > max_angle.degree()){ collapseOK = false; }
    }

    position(v0) = p0;

    // ----- %< -------------------------------------------------------

    // return the result of the collapse simulation
    return collapseOK;
}

float compute_halfedge_priority(pm::vertex_attribute<tg::pos3> const& position, pm::vertex_attribute<Quadric> const& quadrics, pm::halfedge_handle heh)
{
    auto v0(heh.vertex_from());
    auto v1(heh.vertex_to());

    // compute combined quadric that measure distance to all faces incident to either vertex (and those that are already collapsed into them)
    Quadric q = quadrics[v0];
    q += quadrics[v1];

    // evaluate combined quadric for the position of v1 (because the halfedge collapse means v1 will survive)
    return q(position[v1]);
}

void initialize_quadrics(pm::Mesh const& mesh, pm::vertex_attribute<tg::pos3> const& position, pm::face_attribute<tg::vec3> const& normals, pm::vertex_attribute<Quadric>& quadrics)
{
    for (auto vh : mesh.vertices())
    {
        quadrics[vh].clear();

        /*
         * INSERT YOUR OWN CODE BETWEEN THE HORIZONTAL LINES BELOW.
         * DO NOT CHANGE CODE ANYWHERE ELSE.
         *
         * Note that there is another function below this function
         * where you have to insert code as well.
         *
         * Here you should compute each vertex's error quadric as
         * the sum of all fundamental quadrics of their incident
         * faces, respectively. The quadrics are stored as a vertex
         * attribute "quadrics". They are already initialized with
         * the zero-quadric.
         *
         * Quadrics can be constructed by passing the matrix
         * coefficients of the upper right triangle (see
         * QuadricT.hh line 31) or by passing the coefficients
         * of the normal equation defining the plane of the triangle
         * (see Quadric.hh line 43).
         * Quadrics can also be scaled by multiplication with a scalar
         * or combined by adding two Quadrics.
         */

        // ----- %< -------------------------------------------------------

        for (auto face : vh.faces()){
            auto plane_offset = tg::dot(normals(face), position(vh));
            auto distance_to_plane = tg::dot(normals(face), position(vh)) + plane_offset;

            auto const n = normals(face);
            auto const d = distance_to_plane;
            Quadric face_quadric = QuadricT(n.x*n.x, n.x*n.y, n.x*n.z, n.x*d,
                                                     n.y*n.y, n.y*n.z, n.y*d,
                                                              n.z*n.z, n.z*d,
                                                                       d*d);
            quadrics(vh) += face_quadric;
        }

        // ----- %< -------------------------------------------------------
    }
}

void decimate(pm::Mesh& mesh, pm::vertex_attribute<tg::pos3>& position, int num_target_vertices, tg::angle32 max_angle)
{
    auto current_num_vertices(mesh.vertices().count());

    auto priority = mesh.vertices().make_attribute(-1.0f);
    auto collapse_halfedge = mesh.vertices().make_attribute<pm::halfedge_handle>();
    auto quadrics = mesh.vertices().make_attribute(Quadric());
    auto normals = face_normals(position);
    initialize_quadrics(mesh, position, normals, quadrics);

    // compare two vertices based on their priority attribute
    auto cmp = [&](pm::vertex_handle vh1, pm::vertex_handle vh2)
    {
        if (priority[vh1] != priority[vh2])
            return priority[vh1] < priority[vh2];
        return vh1.idx < vh2.idx;
    };

    std::set<pm::vertex_handle, decltype(cmp)> queue(cmp);

    auto enqueue_vertex = [&](pm::vertex_handle vh)
    {
        // find best collapsible halfedge
        pm::halfedge_handle min_hh;

        float prio = -1.0f;
        float min_prio(std::numeric_limits<float>::max());

        for (auto heh : vh.outgoing_halfedges())
        {
            if (!is_collapse_legal(position, normals, heh, max_angle))
                continue;
            prio = compute_halfedge_priority(position, quadrics, heh);
            if (prio != -1.0f && prio < min_prio)
            {
                min_prio = prio;
                min_hh = heh;
            }
        }

        // update queue

        // remove old vertex if present
        if (priority[vh] != -1.0f)
        {
            queue.erase(vh);
            priority[vh] = -1.0;
        }

        // reinsert vertex if there exists a valid collapse
        if (min_hh.is_valid())
        {
            priority[vh] = min_prio;
            collapse_halfedge[vh] = min_hh;
            queue.insert(vh);
        }
    };

    // Build priority queue...
    for (auto vh : mesh.vertices())
        enqueue_vertex(vh);

    // Do the decimation...
    while (current_num_vertices > num_target_vertices && !queue.empty())
    {
        // take first element out of queue
        auto vh = *queue.begin();
        queue.erase(queue.begin());

        /*
         * INSERT YOUR OWN CODE BETWEEN THE HORIZONTAL LINES BELOW.
         * DO NOT CHANGE CODE ANYWHERE ELSE.
         *
         * This is the core of the incremental decimation algorithm.
         * The following steps should be performed:
         *
         * 1.) Perform half-edge collapse operation. Note: Use attribute
         * collapse_halfedge to get the corresponding halfedge
         * that is collapsed for vh). Use function mesh.halfedges().collapse(halfedge_handle)
         * to perform a half-edge collapse. Check if a collapse is legal
         * using the function is_collapse_legal(...) which you modified earlier.
         *
         * 2.) Update the error quadric by adding the quadrics of the
         * collapsing vertices.
         *
         * 3.) Locally update queue by calling enqueue_vertex for every vertex
         * for which the priority might have changed by the latest collapse.
         * enqueue_vertex updates the position of the supplied vertex within the queue as well as
         * storing the best outgoing halfedge in the collapse_halfedge attribute.)
         */

        // ----- %< -------------------------------------------------------


        // ----- %< -------------------------------------------------------
    }

    mesh.compactify();
}


}
