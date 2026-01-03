#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

namespace params = CGAL::parameters;

int main(int argc, char* argv[])
{
  const char* fname = (argc > 1) ? argv[1] : "data/self_intersecting.off";
  std::ifstream input(fname);
  Polyhedron poly;
  if (!input || !(input >> poly)) {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return 1;
  }

  Mesh_domain domain(poly);

  // Rigorous criteria for the official unit test
  Mesh_criteria criteria(params::edge_distance = 0.1,
                         params::edge_min_size = 0.05,
                         params::facet_angle = 25,
                         params::facet_size = 0.1,
                         params::facet_distance = 0.001);

  std::cout << "Starting mesh generation..." << std::endl;
  
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  std::cout << "Finished! Vertices: " << c3t3.triangulation().number_of_vertices() << std::endl;

  return 0;
}
