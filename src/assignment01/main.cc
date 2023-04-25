#include <iostream>
#include <queue>

#include <imgui/imgui.h>

#include <glow-extras/glfw/GlfwContext.hh>
#include <glow-extras/viewer/canvas.hh>
#include <glow-extras/viewer/experimental.hh>

#include <polymesh/Mesh.hh>
#include <polymesh/objects/cone.hh>
#include <polymesh/properties.hh>

#include <typed-geometry/functions/objects/triangle.hh>
#include <typed-geometry/tg.hh>

#include "task.hh"

namespace gp
{
bool in_face(tg::pos2 const& point, pm::face_handle f, pm::vertex_attribute<tg::pos2> const& position)
{
    auto const ps = f.vertices().to_array<3>(position);

    auto const v0 = ps[2] - ps[0];
    auto const v1 = ps[1] - ps[0];
    auto const v2 = point - ps[0];

    // Compute dot products
    auto const dot00 = dot(v0, v0);
    auto const dot01 = dot(v0, v1);
    auto const dot02 = dot(v0, v2);
    auto const dot11 = dot(v1, v1);
    auto const dot12 = dot(v1, v2);

    // Compute barycentric coordinates
    auto const invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    auto const u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    auto const v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    return (u >= 0) && (v >= 0) && (u + v < 1);
}

void project_to_paraboloid(pm::Mesh const& mesh, pm::vertex_attribute<tg::pos2> const& pos2d, pm::vertex_attribute<tg::pos3>& pos3d)
{
    for (auto vh : mesh.vertices())
    {
        auto p = pos2d[vh];
        pos3d[vh] = tg::pos3(p.x, 0.5 + (p.x * p.x + p.y * p.y), p.y);
    }
}

}

int main(int /*argc*/, char** /*args*/)
{
    glow::glfw::GlfwContext ctx;
    tg::rng rng;

    // initialize mesh
    pm::Mesh mesh;
    pm::vertex_attribute<tg::pos2> position(mesh);
    pm::vertex_attribute<tg::pos3> pos3d(mesh);
    pm::vertex_attribute<tg::pos3> paraboloid_pos3d(mesh);
    pm::vertex_attribute<tg::color3> color(mesh);


    // background cone mesh for voronoi diagram
    pm::Mesh cones;
    pm::vertex_attribute<tg::pos3> cone_pos(cones);
    pm::vertex_attribute<tg::color3> cone_color(cones);

    auto const add_cone = [&](tg::pos3 const& apex_position, tg::color3 const& color)
    {
        constexpr int segments = 50;
        pm::objects::add_cone(
            cones,
            [&](pm::vertex_handle v, float x, float y)
            {
                auto const length = 1.0f;
                auto const radius = 1.0f;
                tg::pos3 pos = {radius * tg::cos(x * 2 * tg::pi<float>), radius * sin(x * 2 * tg::pi<float>), -length};
                pos = tg::lerp(pos, tg::pos3(0), y); // apex
                cone_pos[v] = pos + tg::vec3(apex_position);
                cone_color[v] = color;
            },
            segments, true);
    };

    auto const update_paraboloid = [&]() { gp::project_to_paraboloid(mesh, position, paraboloid_pos3d); };
    update_paraboloid();

    auto const reset = [&]()
    {
        rng = tg::rng();
        mesh.clear();
        cones.clear();
        auto const v0 = mesh.vertices().add();
        auto const v1 = mesh.vertices().add();
        auto const v2 = mesh.vertices().add();
        auto const v3 = mesh.vertices().add();

        position[v0] = tg::pos2(-0.5, -0.5);
        position[v1] = tg::pos2(0.5, -0.5);
        position[v2] = tg::pos2(0.5, 0.5);
        position[v3] = tg::pos2(-0.5, 0.5);

        pos3d[v0] = tg::pos3(position[v0]);
        pos3d[v1] = tg::pos3(position[v1]);
        pos3d[v2] = tg::pos3(position[v2]);
        pos3d[v3] = tg::pos3(position[v3]);

        color[v0] = tg::uniform<tg::color3>(rng);
        color[v1] = tg::uniform<tg::color3>(rng);
        color[v2] = tg::uniform<tg::color3>(rng);
        color[v3] = tg::uniform<tg::color3>(rng);

        mesh.faces().add(v0, v1, v2);
        mesh.faces().add(v0, v2, v3);

        // add initial cones
        for (auto const v : mesh.vertices())
            add_cone(tg::pos3(position[v]), color[v]);

        update_paraboloid();
    };

    reset();

    // app state
    bool show_cones = true;
    bool show_triangulation = true;
    bool show_vertices = true;
    bool show_hovered_circumcircle = false;
    bool mouse_locked = false;
    bool show_paraboloid = false;
    bool lock_camera = true;

    std::vector<tg::segment3> circumcircle;

    auto cam_ortho = gv::CameraController::create();
    cam_ortho->setTransform({0, 0, 1}, {0, 0, -1});

    auto cam_free = gv::CameraController::create();
    cam_free->setupMesh(1, {});
    cam_free->setTransform({0, 0, 1}, {0, 0, -1});

    auto cam_para = gv::CameraController::create();
    cam_para->setupMesh(1, {});
    cam_para->setTransform({2, 2, 2}, {0, 0.5, 0});

    auto cone_r = gv::make_renderable(cone_pos);
    cone_r->addAttribute(gv::detail::make_mesh_attribute("aColor", cone_color));
    auto line_r = gv::make_renderable(gv::lines(pos3d).camera_facing());

    auto para_pos_r = gv::make_renderable(paraboloid_pos3d);
    auto para_lines_r = gv::make_renderable(gv::lines(paraboloid_pos3d).line_width_world(0.001));

    int skip = 2; // dirty hack to get rid of warnings during window initialization
    gv::interactive(
        [&](float)
        {
            auto changed = false;

            ImGui::Begin("Options");
            changed |= ImGui::Checkbox("Show Voronoi", &show_cones);
            changed |= ImGui::Checkbox("Show Delaunay", &show_triangulation);
            changed |= ImGui::Checkbox("Show Vertices", &show_vertices);
            changed |= ImGui::Checkbox("Show Circumcircle", &show_hovered_circumcircle);
            changed |= ImGui::Checkbox("Show Paraboloid", &show_paraboloid);
            changed |= ImGui::Checkbox("Lock Camera", &lock_camera);

            if (ImGui::Button("Reset"))
            {
                changed = true;
                reset();
                cone_r = gv::make_renderable(cone_pos);
                cone_r->addAttribute(gv::detail::make_mesh_attribute("aColor", cone_color));
                line_r = gv::make_renderable(gv::lines(pos3d).camera_facing());
            }
            ImGui::End();

            auto const mouse_pos = gv::experimental::interactive_get_mouse_position();

            if (skip > 0)
            {
                --skip;
            }
            else if (ImGui::IsMouseClicked(0) || show_hovered_circumcircle)
            {
                auto const pick = gv::experimental::interactive_get_position(mouse_pos);
                if (pick.has_value() && !ImGui::GetIO().WantCaptureMouse)
                {
                    auto const pick_pos = tg::pos2(pick.value());
                    pm::face_handle picked_face;
                    for (auto const f : mesh.faces())
                    {
                        if (gp::in_face(pick_pos, f, position))
                        {
                            picked_face = f;
                            break;
                        }
                    }

                    if (picked_face.is_valid())
                    {
                        if (show_hovered_circumcircle)
                        {
                            changed = true;
                            auto const ps = picked_face.vertices().to_array<3>(position);
                            auto const cc = tg::circumcircle_of(tg::triangle2{ps[0], ps[1], ps[2]});
                            auto const n_segs = 64;
                            circumcircle.clear();
                            for (auto i = 0; i < n_segs; ++i)
                            {
                                auto const a0 = float(i) * 360_deg / float(n_segs);
                                auto const a1 = float(i + 1) * 360_deg / float(n_segs);
                                auto const s0 = tg::sin(a0);
                                auto const c0 = tg::cos(a0);
                                auto const s1 = tg::sin(a1);
                                auto const c1 = tg::cos(a1);
                                auto const p0 = tg::pos3(cc.center) + tg::vec3(s0, c0, 0.003) * cc.radius;
                                auto const p1 = tg::pos3(cc.center) + tg::vec3(s1, c1, 0.003) * cc.radius;
                                circumcircle.emplace_back(p0, p1);
                            }
                        }

                        if (ImGui::IsMouseClicked(0))
                        {
                            if (position.all([&](tg::pos2 p) { return p != pick_pos; }))
                            {
                                changed = true;
                                auto const v = task::insert_vertex(mesh, position, pick_pos, picked_face);
                                color[v] = tg::uniform<tg::color3>(rng);
                                pos3d[v] = tg::pos3(pick_pos);
                                add_cone(tg::pos3(position[v]), color[v]);
                                cone_r = gv::make_renderable(cone_pos);
                                cone_r->addAttribute(gv::detail::make_mesh_attribute("aColor", cone_color));
                                line_r = gv::make_renderable(gv::lines(pos3d).camera_facing());
                            }
                        }
                    }
                }
            }

            auto g = gv::grid();

            {
                auto v = gv::view();
                if (lock_camera)
                {
                    cam_ortho->setTransform({0, 0, 1}, {0, 0, -1});
                    cam_ortho->enableOrthographicMode({{-0.53f, -0.53f, 0.01f}, {0.53f, 0.53f, 10.0f}});
                    v.configure(cam_ortho);
                }
                else
                {
                    v.configure(cam_free);
                }

                v.configure(gv::no_shadow);
                v.configure(gv::no_ssao);
                v.configure(gv::no_grid);
                v.configure(gv::clear_accumulation(changed));

                if (show_hovered_circumcircle)
                    gv::view(gv::lines(circumcircle).camera_facing(), tg::color3::red, gv::maybe_empty);

                if (mouse_locked)
                {
                    v.configure(gv::no_left_mouse_control);
                    v.configure(gv::no_right_mouse_control);
                }

                if (show_triangulation)
                    gv::view(line_r);

                if (show_cones)
                    gv::view(cone_r, gv::no_shading, gv::no_outline, gv::no_fresnel);

                if (show_vertices)
                {
                    auto c = gv::canvas();
                    for (auto const& p : pos3d)
                        c.add_splat(p + tg::vec3(0, 0, 0.001), tg::vec3{0, 0, 1}).size(0.002).color(tg::color3::magenta);
                }
            }
            if (show_paraboloid)
            {
                if (changed)
                {
                    update_paraboloid();
                    para_pos_r = gv::make_renderable(paraboloid_pos3d);
                    para_lines_r = gv::make_renderable(gv::lines(paraboloid_pos3d).line_width_world(0.001));
                }
                auto v = gv::view(para_pos_r, cam_para, gv::clear_accumulation(changed), gv::no_shadow,
                                  tg::aabb3(tg::pos3(-0.5, 0, -0.5), tg::pos3(0.5, 0.5, 0.5)));
                gv::view(para_lines_r);
            }
        });
}
