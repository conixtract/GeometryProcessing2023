#include <polymesh/Mesh.hh>
#include <typed-geometry/tg-lean.hh>

namespace task
{
tg::dir3 compute_normal(std::vector<pm::vertex_handle> const& vs, pm::vertex_attribute<tg::pos3> const& position);

float compute_mst_weight(pm::vertex_handle v0, pm::vertex_handle v1, pm::vertex_attribute<tg::pos3> const& position, pm::vertex_attribute<tg::dir3> const& normal);
}
