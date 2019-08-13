#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>

using LCC_3 = CGAL::Linear_cell_complex_for_combinatorial_map<2, 3>;
using CST = CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3>;
using Path_on_surface = CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;

LCC_3 lcc;

int main(int argc, char* argv[])
{
  std::cout << "Program shortest_noncontractible_cycle_using_BFS started.\n";
  std::ifstream inp;
  if (argc == 1) inp = std::ifstream("data/3torus.off");
  else inp = std::ifstream(argv[1]);
  if (inp.fail()) {
    std::cout << "Cannot read file. Exiting program\n";
    return EXIT_FAILURE;
  }
  CGAL::load_off(lcc, inp);
  std::cout << "File loaded. Running the main program...\n";
  
  CST                cst(lcc);

  /// Change the value of `root` to test the algorithm at another vertex
  auto root = lcc.darts().begin();
  std::cout << "Finding the shortest noncontractible cycle...\n";
  Path_on_surface cycle = cst.compute_shortest_noncontractible_cycle_with_basepoint(root);
  if (cycle.length() == 0) {
    std::cout << "  Cannot find such cycle. Stop.\n";
    return 0;
  }
  std::cout << "  Number of edges in cycle: " << cycle.length() << std::endl;
  std::cout << "  Root: " << lcc.point_of_vertex_attribute(lcc.vertex_attribute(root)) << std::endl;
}
