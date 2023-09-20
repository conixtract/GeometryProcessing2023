#include <polymesh/Mesh.hh>
#include <typed-geometry/tg-lean.hh>

namespace task
{

pm::edge_attribute<float> compute_weights(pm::Mesh& mesh, pm::vertex_attribute<tg::pos3>& position, bool cotan_weights);

void compute_new_positions(pm::Mesh& mesh, pm::vertex_attribute<tg::pos3>& position, const pm::edge_attribute<float>& edge_weight, const pm::vertex_attribute<bool>& locked, bool simple_laplace, int iterations);

}
