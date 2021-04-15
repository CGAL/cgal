#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef Mesh::Property_map<face_descriptor, std::size_t> Face_proxy_pmap;

namespace VSA = CGAL::Surface_mesh_approximation;

int main()
{
  // reads input mesh
  Mesh mesh;
  std::ifstream file("data/bear.off");
  file >> mesh;

  // face proxy index property map
  Face_proxy_pmap fpxmap =
    mesh.add_property_map<face_descriptor, std::size_t>("f:proxy_id", 0).first;

  // free function interface with named parameters
  VSA::approximate_triangle_mesh(mesh,
  CGAL::parameters::max_number_of_proxies(200). // first stop criterion
    min_error_drop(0.05). // second stop criterion
    number_of_iterations(30). // number of relaxation iterations after seeding
    face_proxy_map(fpxmap)); // output face-proxy map

  // iterates over faces and outputs segment id to console
  for(face_descriptor f : faces(mesh))
    std::cout << fpxmap[f] << std::endl;

  return EXIT_SUCCESS;
}
