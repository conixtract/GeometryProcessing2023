#include <chrono>
#include <fstream>
#include <iostream>
#include <queue>

#include <imgui/imgui.h>
#include <glow-extras/glfw/GlfwContext.hh>
#include <glow-extras/viewer/view.hh>
#include <polymesh/Mesh.hh>
#include <polymesh/algorithms/normalize.hh>
#include <polymesh/formats.hh>
#include <polymesh/properties.hh>
#include <typed-geometry/tg.hh>

std::vector<tg::segment3> get_arrow(pm::vertex_attribute<tg::pos3> const& position, pm::halfedge_handle heh)
{
    auto p0 = position(heh.vertex_from());
    auto p1 = position(heh.vertex_to());
    auto p2 = position(heh.next().vertex_to());
    auto shaft_start = (0.7f * p0 + 0.15f * p1 + 0.15f * p2) / 1.0;
    auto shaft_end = (0.15f * p0 + 0.7f * p1 + 0.15f * p2) / 1.0;
    auto dir = 0.2f * (shaft_end - shaft_start);
    auto n = tg::normalize(tg::cross(p1 - p0, p2 - p0));
    auto dir_90 = tg::cross(n, dir);
    auto point_start = shaft_end - dir + 0.5 * dir_90;
    auto point_end = shaft_end;

    std::vector<tg::segment3> arrow{tg::segment3(shaft_start, shaft_end), tg::segment3(point_start, point_end)};

    return arrow;
}

auto view(pm::vertex_attribute<tg::pos3> const& position, std::vector<pm::vertex_handle> const& vertices, tg::color3 const& color)
{
    std::vector<tg::pos3> points;
    for (auto vh : vertices)
        points.emplace_back(position[vh]);
    return gv::view(gv::points(points).point_size_world(0.01f), color, gv::maybe_empty);
}

auto view(pm::vertex_attribute<tg::pos3> const& position,
          pm::halfedge_handle heh,
          pm::vertex_handle center,
          std::vector<pm::vertex_handle> const& reached_vertices,
          std::vector<pm::vertex_handle> const& to_be_reached)
{
    auto v = gv::view();
    //  current halfedge
    if (heh.is_valid())
        gv::view(gv::lines(get_arrow(position, heh)).line_width_world(0.007f), tg::color3::blue);

    // vertices to be reached
    view(position, to_be_reached, tg::color3::red);

    // vertices already reached
    view(position, reached_vertices, tg::color3::green);

    // center vertex
    if (center.is_valid())
        gv::view(tg::sphere3(position[center], 0.01f), tg::color3::magenta);

    return v;
}

bool update_targets(std::vector<pm::vertex_handle>& to_be_reached, std::vector<pm::vertex_handle>& reached, pm::halfedge_handle heh)
{
    bool changed = false;
    auto it = std::find(to_be_reached.begin(), to_be_reached.end(), heh.vertex_to());
    if (it != to_be_reached.end())
    {
        reached.push_back(*it);
        to_be_reached.erase(it);
        changed = true;
    }
    return changed;
}

void game(pm::vertex_attribute<tg::pos3> const& position, pm::halfedge_handle selected_heh, pm::vertex_handle center_vertex, std::vector<pm::vertex_handle>& to_be_reached)
{
    std::vector<pm::vertex_handle> reached;
    auto rm = gv::make_renderable(position);
    auto l = gv::make_renderable(gv::lines(position).line_width_world(0.005f));
    gv::interactive([&](float) {
        bool changed = false;
        ImGui::Begin("Navigation");
        if (ImGui::Button("Next"))
        {
            selected_heh = selected_heh.next();
            changed = true;
        }
        if (ImGui::Button("Opposite"))
        {
            selected_heh = selected_heh.opposite();
            changed = true;
        }

        ImGui::End();

        changed |= update_targets(to_be_reached, reached, selected_heh);
        auto v = gv::view(rm, gv::clear_accumulation(changed));
        gv::view(l);
        view(position, selected_heh, center_vertex, reached, to_be_reached);

        if (to_be_reached.empty())
            gv::close_viewer();
    });
}

void game1(pm::Mesh const& mesh, pm::vertex_attribute<tg::pos3> const& position)
{
    auto circulator_center = mesh.faces()[1561];
    auto to_be_reached = circulator_center.vertices().to_vector();
    auto selected_heh = circulator_center.halfedges().first();

    game(position, selected_heh, pm::vertex_handle(), to_be_reached);
}

void game2(pm::Mesh const& mesh, pm::vertex_attribute<tg::pos3> const& position)
{
    auto circulator_center = mesh.vertices()[1361];
    auto to_be_reached = circulator_center.adjacent_vertices().to_vector();
    auto selected_heh = circulator_center.incoming_halfedges().first();

    game(position, selected_heh, circulator_center, to_be_reached);
}

void game3(pm::Mesh const& mesh, pm::vertex_attribute<tg::pos3> const& position)
{
    std::vector<pm::vertex_handle> to_be_reached{mesh.vertices()[804], mesh.vertices()[1775]};
    auto selected_heh = to_be_reached.front().incoming_halfedges().first();

    game(position, selected_heh, pm::vertex_handle(), to_be_reached);
}

int main(int /*argc*/, char** /*args*/)
{
    glow::glfw::GlfwContext ctx;

    // initialize mesh
    pm::Mesh mesh;
    auto position = mesh.vertices().make_attribute<tg::pos3>();

    // load the mesh
    auto const search_paths = std::array{
        "../../../src/assignment00/data/bunny_open.obj", //
        "../../src/assignment00/data/bunny_open.obj",    //
        "../src/assignment00/data/bunny_open.obj"        //
    };

    // check if file exists.
    // note: don't use std::filesystem::exists because macOS supports just it since 10.15 (October 2019)
    auto const exists = [](std::string const& path) {
        std::ifstream f(path.c_str());
        return f.good();
    };

    for (auto const p : search_paths)
    {
        if (exists(p))
            pm::load(p, mesh, position);
    }
    if (mesh.vertices().empty())
    {
        std::cerr << "Failed to read input mesh. Please set your working directory to ./bin" << std::endl;
        return -1;
    }

    pm::normalize(position);

    auto begin = std::chrono::high_resolution_clock::now();

    game1(mesh, position);
    auto end1 = std::chrono::high_resolution_clock::now();

    std::cout << "Visiting all vertices of a face took you " << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - begin).count() << "ms." << std::endl;

    game2(mesh, position);
    auto end2 = std::chrono::high_resolution_clock::now();

    std::cout << "Visiting the one ring of a vertex took you " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end1).count() << "ms."
              << std::endl;

    game3(mesh, position);
    auto end3 = std::chrono::high_resolution_clock::now();

    std::cout << "Navigating across the mesh took you " << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - end2).count() << "ms." << std::endl;

    std::cout << "Your total time is " << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - begin).count() << "ms!" << std::endl;
}
