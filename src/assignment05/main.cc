#include <fstream>

#include <imgui/imgui.h>
#include <glow-extras/glfw/GlfwContext.hh>
#include <glow-extras/viewer/view.hh>
#include <polymesh/Mesh.hh>
#include <polymesh/algorithms/normalize.hh>
#include <polymesh/formats.hh>
#include <typed-geometry/tg-lean.hh>

#include "task.hh"

namespace
{
bool exists(std::string const& filepath)
{
    std::ifstream infile(filepath);
    return infile.good();
}
}

int main(int /*argc*/, char** /*args*/)
{
    // used for rendering
    glow::glfw::GlfwContext ctx;

    pm::Mesh mesh;
    // the vertex positions
    pm::vertex_attribute<tg::pos3> position(mesh);
    // 3d position from 2d texture coordinates
    pm::vertex_attribute<tg::pos3> position_from_tex_coord(mesh);
    // texture coordinates
    pm::vertex_attribute<tg::pos2> texture_coordinates(mesh);
    // edge weights (uniform / cotangent)
    pm::edge_attribute<float> edge_weight(mesh);
    // parameter domain as color attribute
    pm::vertex_attribute<tg::color3> paramter_color(mesh);
    // use world x and y as texture coordinate
    pm::vertex_attribute<tg::pos2> texture_from_world(mesh);

    // renderables
    decltype(gv::make_renderable(position)) position_r;
    decltype(gv::make_renderable(gv::lines(position))) wireframe_r;
    decltype(gv::make_renderable(position)) position_circle_r;
    decltype(gv::make_renderable(gv::lines(position))) wireframe_circle_r;

    // folders to search for data
    std::array<char const*, 3> folders = {"../../../src/assignment05/data/", "../../src/assignment05/data/", "../src/assignment05/data/"};
    std::array<const char*, 3> mesh_filenames = {"bunny_open.off", "Circle.off", "hemisphere.off"};
    std::string texture_filename = "gp-checker.png";

    std::string folder;
    for (auto const f : folders)
        if (exists(f + texture_filename))
        {
            folder = f;
            break;
        }

    if (folder.empty())
    {
        std::cerr << "Wrong working directory! Make sure your current working directory is ./bin!" << std::endl;
        exit(1);
    }

    // checkerboard texture
    glow::SharedTexture2D tex = glow::Texture2D::createFromFile(folder + texture_filename, glow::ColorSpace::sRGB);

    // gui state
    bool show_parameter_domain = false;
    bool use_world_position_as_tex_coord = false;
    gp::weight_type selected_weight_type = gp::weight_type::uniform;
    const char* weight_type_name[] = {"uniform", "cotangent"};
    int selected_mesh_idx = 0;
    int iterations = 10;

    // update renderables whenever they are changed
    auto update_renderables = [&]() {
        position_r = gv::make_renderable(position);

        auto const& tex_coord = use_world_position_as_tex_coord ? texture_from_world : texture_coordinates;

        for (auto const v : mesh.vertices())
            position_from_tex_coord[v] = tg::pos3(texture_coordinates[v]);

        position_circle_r = gv::make_renderable(position_from_tex_coord);
        wireframe_circle_r = gv::make_renderable(gv::lines(position_from_tex_coord).line_width_world(0.002));

        if (show_parameter_domain)
        {
            for (auto const v : mesh.vertices())
                paramter_color[v] = {tex_coord[v].x, tex_coord[v].y, 0};
            gv::configure(*position_r, paramter_color);
            gv::configure(*position_circle_r, paramter_color);
        }
        else
        {
            gv::configure(*position_r, gv::textured(tex_coord, tex));
            gv::configure(*position_circle_r, gv::textured(tex_coord, tex));
        }

        wireframe_r = gv::make_renderable(gv::lines(position).line_width_world(0.002));
    };

    auto const load_mesh = [&](std::string const& filename) {
        mesh.clear();
        pm::load(folder + filename, mesh, position);
        pm::normalize(position);
        for (auto& p : position)
        {
            p *= 0.5f;
            p += tg::vec3(0.5, 0.5, 0);
        }
        task::init_texture_coordinates(position, texture_coordinates);
        for (auto const v : mesh.vertices())
            texture_from_world[v] = {position[v].x, position[v].y};
        task::compute_weights(selected_weight_type, position, edge_weight);
    };

    load_mesh(mesh_filenames[selected_mesh_idx]);

    update_renderables();

    // viewer loop
    gv::interactive([&](auto) {
        auto changed = false;
        ImGui::Begin("Parametrization");
        if (ImGui::Combo("Mesh", &selected_mesh_idx, mesh_filenames.data(), 3))
        {
            load_mesh(mesh_filenames[selected_mesh_idx]);
            changed |= true;
        }
        if (ImGui::Combo("Weight Type", reinterpret_cast<int*>(&selected_weight_type), weight_type_name, 2))
        {
            task::compute_weights(selected_weight_type, position, edge_weight);
        }
        ImGui::InputInt("Iterations", &iterations);
        if (ImGui::Button("Restart"))
        {
            task::init_texture_coordinates(position, texture_coordinates);
            changed |= true;
        }
        if (ImGui::Button("Smooth"))
        {
            task::smooth_texcoords(mesh, iterations, edge_weight, texture_coordinates);
            changed |= true;
        }
        if (ImGui::Button("Parametrize"))
        {
            task::direct_solve(position, edge_weight, selected_weight_type, texture_coordinates);
            changed |= true;
        }

        changed |= ImGui::Checkbox("Show Parameter Domain", &show_parameter_domain);
        changed |= ImGui::Checkbox("Use world position for texture coordinates", &use_world_position_as_tex_coord);

        ImGui::End();

        if (changed)
            update_renderables();

        auto v = gv::grid();
        {
            auto v = gv::view(position_r, gv::clear_accumulation(changed));
            gv::view(wireframe_r);
        }
        {
            auto v = gv::view(position_circle_r, gv::clear_accumulation(changed));
            gv::view(wireframe_circle_r);
        }
    });
}
