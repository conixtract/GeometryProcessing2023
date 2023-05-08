#include <algorithm>
#include <iostream>
#include <queue>

#include <imgui/imgui.h>
#include <glow-extras/glfw/GlfwContext.hh>
#include <glow-extras/viewer/view.hh>
#include <polymesh/Mesh.hh>
#include <polymesh/algorithms/normalize.hh>
#include <polymesh/formats.hh>
#include <typed-geometry/tg.hh>

#include "task.hh"

namespace gp
{
std::vector<pm::vertex_handle> k_nearest_neighbors(pm::vertex_handle v, pm::vertex_attribute<tg::pos3> const& position, int k)
{
    auto const safe_k = tg::min(k, position.mesh().vertices().size());

    // all vertices
    auto neighbors = v.mesh->vertices().to_vector();

    // remove v from its neighbors
    neighbors.erase(std::remove(neighbors.begin(), neighbors.end(), v), neighbors.end());

    // sort by distance from v
    std::partial_sort(neighbors.begin(), neighbors.begin() + safe_k, neighbors.end(), [&](pm::vertex_handle v0, pm::vertex_handle v1) {
        return tg::distance_sqr(position[v], position[v0]) < tg::distance_sqr(position[v], position[v1]);
    });

    // keep only k
    neighbors.resize(safe_k);

    return neighbors;
}

void propagate_orientation(pm::vertex_handle v_from, pm::vertex_handle v_to, pm::vertex_attribute<tg::dir3>& normal)
{
    // flip the normal of v_to, if it points the opposite way of v_from's normal
    if (dot(normal[v_from], normal[v_to]) < 0.0f)
        normal[v_to] = -normal[v_to];
}

void propagate_orientation(pm::vertex_attribute<tg::pos3> const& position, int k, pm::vertex_handle seed, pm::vertex_attribute<tg::dir3>& normal, std::vector<std::pair<pm::vertex_handle,pm::vertex_handle>>& spanning_tree_edges)
{
    // candidate edges in the minimum spanning tree
    struct propagation_candidate
    {
        pm::vertex_handle v_from;
        pm::vertex_handle v_to;
        float edge_weight;
        // std::priority_queue takes the largest element by default.
        // Since we want the edge with the minimal weight, we need to invert the comparison!
        bool operator<(propagation_candidate const& rhs) const { return edge_weight > rhs.edge_weight; }
    };

    auto const& mesh = position.mesh();

    // visited flag. Each vertex should be visited only once
    auto visited = mesh.vertices().make_attribute(false);

    // initialize nearest neighbor graph
    pm::vertex_attribute<std::vector<pm::vertex_handle>> neighbors(mesh);
    for (auto const v : position.mesh().vertices())
        neighbors[v] = k_nearest_neighbors(v, position, k);

    // priority queue to build the minimum spanning tree
    std::priority_queue<propagation_candidate> q;

    // initialize with a (given) random vertex
    q.push({seed, seed, 0.0f});

    spanning_tree_edges.clear();
    while (!q.empty())
    {
        // get the top element
        auto const current = q.top();
        q.pop();

        // skip vertices that were already visited
        if (visited[current.v_to])
          continue;

        // set vertex visited
        visited[current.v_to] = true;

        // propagate the normal orientation
        propagate_orientation(current.v_from, current.v_to, normal);
        spanning_tree_edges.push_back({current.v_from, current.v_to});

        // add neighbors to priority queue
        for (auto const& neighbor : neighbors[current.v_to])
        {
            // only visit each vertex once
            if (!visited[neighbor])
            {
                q.push({current.v_to, neighbor, task::compute_mst_weight(current.v_to, neighbor, position, normal)});
            }
        }
    }
}


void estimate_normals(pm::Mesh const& mesh, pm::vertex_attribute<tg::pos3> const& position, int k, pm::vertex_attribute<tg::dir3>& normal)
{
    auto count = 0;
    for (auto const v : mesh.vertices())
    {
        auto const neighbors = k_nearest_neighbors(v, position, k);
        normal[v] = task::compute_normal(neighbors, position);

        // report progress
        std::cout << "vertex " << ++count << " of " << mesh.vertices().size() << std::endl;
    }
}
}

int main(int /*argc*/, char** /*args*/)
{
    // used for rendering
    glow::glfw::GlfwContext ctx;

    pm::Mesh mesh;
    // the vertex positions
    pm::vertex_attribute<tg::pos3> position(mesh);
    // the vertex normals
    pm::vertex_attribute<tg::dir3> normal(mesh);
    // the k in k-nearest-neighbors
    int const k = 10;

    // where to start spanning tree
    pm::vertex_handle seed;

    tg::aabb3 aabb;

    std::vector<std::pair<pm::vertex_handle,pm::vertex_handle>> spanning_tree_edges;

    // folders to search for data
    const char* folders[] = {"../../../src/assignment02/data/", //
                             "../../src/assignment02/data/",    //
                             "../src/assignment02/data/"};

    char const* filenames[] = {"sphere.off", "tetra_thing.off"};

    auto const load = [&](std::string const& filename) {
        mesh.clear();
        if (!pm::load(folders[0] + filename, mesh, position) && //
            !pm::load(folders[1] + filename, mesh, position) && //
            !pm::load(folders[2] + filename, mesh, position))
        {
            std::cerr << "Failed to load file. Make sure your current working directory is gp2020assignments/bin!" << std::endl;
            exit(1);
        }
        pm::normalize(position);
        aabb = tg::aabb_of(position);
        seed = mesh.vertices().first();
    };
    load(filenames[0]);

    // gui state
    int current_item = 0;
    bool normals_computed = false;
    bool spanning_tree_computed = false;
    bool show_normals = true;
    bool show_spanning_tree = false;
    std::vector<tg::segment3> normal_segments;
    std::vector<tg::segment3> spanning_tree_segments;

    auto const compute_normal_segments = [&]() {
        normal_segments.clear();
        for (auto const v : mesh.vertices())
            normal_segments.push_back({position[v], position[v] + normal[v] * 0.07f});
    };
    auto const compute_spanning_tree_segments = [&]() {
        spanning_tree_segments.clear();
        for (auto const& pair : spanning_tree_edges)
            spanning_tree_segments.push_back({position[pair.first], position[pair.second]});
    };

    auto point_r = gv::make_renderable(gv::points(position).point_size_world(0.01));
    auto normal_r = gv::make_renderable(gv::lines(normal_segments).line_width_world(0.005));
    auto spanning_tree_r = gv::make_renderable(gv::lines(normal_segments).line_width_world(0.005));
    gv::configure(*normal_r, tg::color3::red);
    gv::configure(*spanning_tree_r, tg::color3::green);

    gv::interactive([&](auto) {
        ImGui::Begin("Normal Estimation");
        auto data_changed = ImGui::Combo("Data", &current_item, filenames, 2);
        auto normals_changed = false;
        if (ImGui::Button("Compute Normals"))
        {
            std::cout << "Estimating normals" << std::endl;
            gp::estimate_normals(mesh, position, k, normal);
            compute_normal_segments();
            normals_computed = true;
            normals_changed = true;
            normal_r = gv::make_renderable(gv::lines(normal_segments).line_width_world(0.005));
            gv::configure(*normal_r, tg::color3::red);
        }
        if (ImGui::Button("Propagate Normals"))
        {
            if (!normals_computed)
            {
                std::cout << "Estimating normals" << std::endl;
                gp::estimate_normals(mesh, position, k, normal);
                compute_normal_segments();
                normals_computed = true;
            }

            // start with the vertex with the smallest value in x direction
            for (auto const v : mesh.vertices())
            {
                if (position[v].x < position[seed].x)
                    seed = v;
            }
            // flip seed normal if necessary
            if (normal[seed].x > 0.0f)
                normal[seed] = -normal[seed];

            gp::propagate_orientation(position, k, seed, normal, spanning_tree_edges);

            compute_normal_segments();
            normals_changed = true;
            normal_r = gv::make_renderable(gv::lines(normal_segments).line_width_world(0.005));
            gv::configure(*normal_r, tg::color3::red);

            compute_spanning_tree_segments();
            spanning_tree_r = gv::make_renderable(gv::lines(spanning_tree_segments).line_width_world(0.015));
            gv::configure(*spanning_tree_r, tg::color3::green);
            spanning_tree_computed = true;
        }
        bool something_changed = ImGui::Checkbox("Show normals", &show_normals);
        something_changed |= ImGui::Checkbox("Show Spanning Tree", &show_spanning_tree);
        ImGui::End();

        if (data_changed)
        {
            load(filenames[current_item]);
            normals_computed = false;
            spanning_tree_computed = false;
            point_r = gv::make_renderable(gv::points(position).point_size_world(0.01));
        }

        auto v = gv::view(point_r, aabb, gv::clear_accumulation(normals_changed || data_changed || something_changed));
        if (normals_computed && show_normals)
            gv::view(normal_r);

        if (spanning_tree_computed && show_spanning_tree)
        {
            gv::view(tg::sphere3(position[seed], 0.02f), tg::color3::magenta);
            gv::view(spanning_tree_r);
        }
    });
}
