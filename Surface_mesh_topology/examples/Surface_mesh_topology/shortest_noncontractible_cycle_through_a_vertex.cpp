#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/squared_distance_3.h>

using LCC_3             = CGAL::Linear_cell_complex_for_generalized_map<2, 3>;
using Dart_const_handle = LCC_3::Dart_const_handle;
using Path_on_surface   = CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;
using CST               = CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3>;

struct Weight_functor {
  Weight_functor(const LCC_3& lcc) : m_lcc(lcc) { }
  using Weight_t = double;
  Weight_t operator()(Dart_const_handle dh) const {
    auto x = m_lcc.point_of_vertex_attribute(m_lcc.vertex_attribute(dh));
    auto y = m_lcc.point_of_vertex_attribute(m_lcc.vertex_attribute(m_lcc.template alpha<0>(dh)));
    return CGAL::sqrt(CGAL::squared_distance(x, y));
  }
private:
  const LCC_3& m_lcc;
};

LCC_3 lcc;

int main(int argc, char* argv[])
{
  std::cout << "Program shortest_noncontractible_cycle_through_a_vertex started.\n";
  std::ifstream inp;
  if (argc == 1) inp = std::ifstream("data/3torus.off");
  else inp = std::ifstream(argv[1]);
  if (inp.fail()) {
    std::cout << "Cannot load file. Exiting program...\n";
    return EXIT_FAILURE;
  }
  CGAL::load_off(lcc, inp);
  std::cout << "File loaded. Running the main program...\n";

  CST cst(lcc);
  Weight_functor wf(lcc);

  /// Change the value of `root` to test the algorithm at another vertex
  auto root = lcc.darts().begin();
  std::cout << "Finding the shortest noncontractible cycle...\n";
  Path_on_surface cycle = cst.compute_shortest_noncontractible_cycle_with_basepoint(root, wf);
  if (cycle.length() == 0) {
    std::cout << "  Cannot find such cycle. Stop.\n";
    return 0;
  }
  double cycle_length = 0;
  for (int i = 0; i < cycle.length(); ++i)
    cycle_length += wf(cycle[i]);
  std::cout << "  Number of edges in cycle: " << cycle.length() << std::endl;
  std::cout << "  Cycle length: " << cycle_length << std::endl;
  std::cout << "  Root: " << lcc.point_of_vertex_attribute(lcc.vertex_attribute(root)) << std::endl;
}
