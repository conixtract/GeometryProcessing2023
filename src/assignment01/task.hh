#include <polymesh/Mesh.hh>

#include <typed-geometry/types/pos.hh>

namespace task
{
bool is_delaunay(pm::edge_handle edge, pm::vertex_attribute<tg::pos2> const& position);

pm::vertex_index insert_vertex(pm::Mesh& mesh, pm::vertex_attribute<tg::pos2>& position, tg::pos2 const& vertex_position, pm::face_handle face);
}
