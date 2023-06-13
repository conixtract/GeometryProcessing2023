#include <polymesh/Mesh.hh>
#include <typed-geometry/tg-lean.hh>

namespace task
{
bool is_collapse_legal(pm::vertex_attribute<tg::pos3>& position, pm::face_attribute<tg::vec3> const& normal, pm::halfedge_handle heh, tg::angle32 max_angle);

void decimate(pm::Mesh& mesh, pm::vertex_attribute<tg::pos3>& position, int target_num_vertices, tg::angle32 max_angle);

}
