#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <CGAL/Shortest_noncontractible_cycle.h>
#include <CGAL/IO/Color.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/squared_distance_3.h>
#define CGAL_USE_BASIC_VIEWER 1

using LCC_3             = CGAL::Linear_cell_complex_for_generalized_map<2, 3>;
using Dart_const_handle = LCC_3::Dart_const_handle;

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

using SNC = CGAL::Surface_mesh_topology::Shortest_noncontractible_cycle<LCC_3, Weight_functor>;

LCC_3 lcc;

int main(int argc, char* argv[])
{
  std::cout << "Program shortest_noncontractible_cycle_through_a_vertex started.\n";
  std::ifstream inp;
  if (argc == 1) inp = std::ifstream("../../examples/Surface_mesh_topology/data/3torus.off");
  else inp = std::ifstream(argv[1]);
  if (inp.fail()) {
    std::cout << "Cannot load file. Exiting program...\n";
    return EXIT_FAILURE;
  }
  CGAL::load_off(lcc, inp);
  std::cout << "File loaded. Running the main program...\n";

  Weight_functor     wf(lcc);
  SNC                snc(lcc, wf);
  SNC::Path          cycle;
  SNC::Distance_type cycle_length;

  /// Change the value of `root` to test the algorithm at another vertex
  auto root = lcc.darts().begin();
  std::cout << "Finding the shortest noncontractible cycle...\n";
  snc.find_cycle(root, cycle, &cycle_length);
  if (cycle.size() == 0) {
    std::cout << "  Cannot find such cycle. Stop.\n";
    return 0;
  }
  std::cout << "  Number of edges in cycle: " << cycle.size() << std::endl;
  std::cout << "  Cycle length: " << cycle_length << std::endl;
  std::cout << "  Root: " << lcc.point_of_vertex_attribute(lcc.vertex_attribute(root)) << std::endl;
}
