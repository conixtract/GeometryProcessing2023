#include <polymesh/Mesh.hh>
#include <typed-geometry/tg-lean.hh>

namespace gp
{
enum weight_type
{
    uniform = 0,
    cotangent
};
}

namespace task
{
void init_texture_coordinates(pm::vertex_attribute<tg::pos3> const& position, pm::vertex_attribute<tg::pos2>& texture_coordintate);
void compute_weights(gp::weight_type type, pm::vertex_attribute<tg::pos3> const& position, pm::edge_attribute<float>& edge_weight);
void direct_solve(pm::vertex_attribute<tg::pos3> const& position, pm::edge_attribute<float>& edge_weight, gp::weight_type type, pm::vertex_attribute<tg::pos2>& texture_coordinate);
void smooth_texcoords(pm::Mesh const& m, int iterations, pm::edge_attribute<float> const& weight, pm::vertex_attribute<tg::pos2>& texture_coordinate);
}
