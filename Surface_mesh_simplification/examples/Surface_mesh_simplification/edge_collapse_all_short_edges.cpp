#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>                          Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const char* filename = (argc > 1) ? argv[1] : "data/cube-meshed.off";
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

  const double threshold = (argc > 2) ? std::atof(argv[2]) : 0.2;

  std::cout << "Collapsing edges with length larger than " << threshold
            << " of mesh: " << filename << "..." << std::endl;
  int r = SMS::edge_collapse(surface_mesh,
                             CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<double>(threshold),
                             CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, surface_mesh))
                                              .halfedge_index_map(get(CGAL::halfedge_external_index  ,surface_mesh))
                                              .get_cost(SMS::Edge_length_cost <Surface_mesh>())
                                              .get_placement(SMS::Midpoint_placement<Surface_mesh>()));

  std::cout << "\nFinished!\n" << r << " edges removed.\n"
            << (surface_mesh.size_of_halfedges()/2) << " final edges.\n";

  std::ofstream os(argc > 3 ? argv[3] : "out.off");
  os.precision(17);
  os << surface_mesh;

  return EXIT_SUCCESS;
}
