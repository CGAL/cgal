#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const char* filename = (argc > 1) ? argv[1] : "data/fold.off";
  std::ifstream is(filename);
  if(!is || !(is >> surface_mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // This is a stop predicate (defines when the algorithm terminates).
  // In this example, the simplification stops when the number of undirected edges
  // left in the surface mesh drops below the specified number
  const std::size_t stop_n = (argc > 2) ? std::stoi(argv[2]) : num_halfedges(surface_mesh)/2 - 1;
  SMS::Count_stop_predicate<Surface_mesh> stop(stop_n);

  typedef SMS::Bounded_normal_change_placement<SMS::LindstromTurk_placement<Surface_mesh> > Placement;

  // This the actual call to the simplification algorithm.
  // The surface mesh and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface mesh lack an "id()" field.
  std::cout << "Collapsing edges of mesh: " << filename << ", aiming for " << stop_n << " final edges..." << std::endl;
  SMS::edge_collapse(surface_mesh, stop,
                     CGAL::parameters::get_cost(SMS::LindstromTurk_cost<Surface_mesh>())
                                      .get_placement(Placement()));

  std::ofstream os((argc > 3) ? argv[3] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
