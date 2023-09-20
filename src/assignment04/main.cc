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
}

int main(int /*argc*/, char** /*args*/)
{
    // used for rendering
    glow::glfw::GlfwContext ctx;

    pm::Mesh mesh;
    // the vertex positions
    pm::vertex_attribute<tg::pos3> position(mesh);
    decltype(gv::make_renderable(position)) position_r;
    decltype(gv::make_renderable(gv::lines(position))) wireframe_r;

    auto update_renderables = [&]()
    {
      position_r = gv::make_renderable(position);
      wireframe_r = gv::make_renderable(gv::lines(position).line_width_world(0.002));
    };

    tg::aabb3 aabb;

    // folders to search for data
    const char* folders[] = {"../../../src/assignment04/data/", //
                             "../../src/assignment04/data/",    //
                             "../src/assignment04/data/"};

    std::string filename = "bunny.off";


    // gui state
    int target_num_vertices = 0;
    int max_angle_int = 45;


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
        update_renderables();
        target_num_vertices = static_cast<int>(std::pow(2, 13));
    };
    load(filename);

    gv::interactive([&](auto) {
        ImGui::Begin("Decimation");
        ImGui::InputInt("target #vertices", &target_num_vertices);
        ImGui::InputInt("Max normal deviation", &max_angle_int);
        bool changed = false;
        if (ImGui::Button("Decimate"))
        {
            std::cout << "Decimating" << std::endl;
            task::decimate(mesh, position, target_num_vertices, tg::degree(max_angle_int));
            target_num_vertices /= 4;

            update_renderables();
            changed = true;
        }
        if (ImGui::Button("Reset"))
        {
          load(filename);
          changed = true;
        }
        ImGui::End();

        auto v = gv::view(position_r, gv::clear_accumulation(changed), aabb, "Bunny, #vertices = " + std::to_string(mesh.vertices().count()));
        gv::view(wireframe_r);

    });
}
