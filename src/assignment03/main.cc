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
struct Versions
{
    Versions(std::string name, bool cotan, bool bilaplacian, pm::Mesh& mesh)
      : name(std::move(name)), positions(mesh), weights(mesh), locked(mesh), cotan(cotan), bilaplacian(bilaplacian)
    {
    }

    std::string name;
    pm::vertex_attribute<tg::pos3> positions;
    pm::edge_attribute<float> weights;
    pm::vertex_attribute<bool> locked;
    bool cotan;
    bool bilaplacian;
    decltype(gv::make_renderable(positions)) positions_r;
    decltype(gv::make_renderable(gv::lines(positions))) wireframe_r;
};

}

int main(int /*argc*/, char** /*args*/)
{
    // used for rendering
    glow::glfw::GlfwContext ctx;

    pm::Mesh mesh;
    // the vertex positions
    pm::vertex_attribute<tg::pos3> position_initial(mesh);

    tg::aabb3 aabb;

    std::vector<gp::Versions> versions;
    versions.push_back({"Uniformly Weighted Laplacian Smoothing", false, false, mesh});
    versions.push_back({"Cotan Weighted Laplacian Smoothing", true, false, mesh});
    versions.push_back({"Uniformly Weighted Bi-Laplacian Smoothing", false, true, mesh});
    versions.push_back({"Cotan Weighted Bi-Laplacian Smoothing", true, true, mesh});

    // folders to search for data
    const char* folders[] = {"../../../src/assignment03/data/", //
                             "../../src/assignment03/data/",    //
                             "../src/assignment03/data/"};

    char const* filenames[] = {"blend.off", "blend2.off", "fan.off", "fandisk_noise.off"};


    // gui state
    int current_model = 0;
    int num_iterations = 1;
    bool show_wireframe = false;

    auto const lock_vertices = [&](gp::Versions& v) {
        // always lock boundary
        for (auto vh : mesh.vertices())
            v.locked(vh) = vh.is_boundary();
        if (v.bilaplacian) // for bi-laplacian smoothing also lock one-ring of boundary
            for (auto vh : mesh.vertices())
                v.locked(vh) |= vh.adjacent_vertices().any([&](auto vh) { return vh.is_boundary(); });
    };

    auto const load = [&](std::string const& filename) {
        mesh.clear();
        if (!pm::load(folders[0] + filename, mesh, position_initial) && //
            !pm::load(folders[1] + filename, mesh, position_initial) && //
            !pm::load(folders[2] + filename, mesh, position_initial))
        {
            std::cerr << "Failed to load file. Make sure your current working directory is gp2020assignments/bin!" << std::endl;
            exit(1);
        }
        pm::normalize(position_initial);
        for (auto vh : mesh.vertices())
        {
            auto& p = position_initial(vh);
            p = tg::pos3(-p.y, -p.z, p.x);
        }
        aabb = tg::aabb_of(position_initial);
        for (auto& ver : versions)
        {
            ver.positions = position_initial;
            ver.positions_r = gv::make_renderable(ver.positions);
            ver.wireframe_r = gv::make_renderable(gv::lines(ver.positions).line_width_world(0.005));
            ver.weights = task::compute_weights(mesh, position_initial, ver.cotan);
            lock_vertices(ver);
        }
    };
    load(filenames[0]);

    gv::interactive([&](auto) {
        ImGui::Begin("Smoothing");
        auto data_changed = ImGui::Combo("Data", &current_model, filenames, 4);
        ImGui::InputInt("# iterations", &num_iterations);
        auto wireframe_changed = ImGui::Checkbox("Show Wireframe", &show_wireframe);
        auto positions_changed = false;
        if (ImGui::Button("Smooth"))
        {
            std::cout << "Smoothing" << std::endl;
            for (auto& ver : versions)
            {
                auto edge_weight = task::compute_weights(mesh, position_initial, ver.cotan);
                task::compute_new_positions(mesh, ver.positions, edge_weight, ver.locked, !ver.bilaplacian, num_iterations);
                ver.positions_r = gv::make_renderable(ver.positions);
                ver.wireframe_r = gv::make_renderable(gv::lines(ver.positions).line_width_world(0.005));
            }
            positions_changed = true;
        }
        data_changed |= (ImGui::Button("Reset"));
        ImGui::End();

        if (data_changed)
        {
            load(filenames[current_model]);
        }

        auto g = gv::grid();

        for (const auto& ver : versions)
        {
            auto v = gv::view(ver.positions_r, gv::clear_accumulation(positions_changed || data_changed || wireframe_changed), aabb, ver.name);
            if (show_wireframe)
                gv::view(ver.wireframe_r);
        }
    });
}
