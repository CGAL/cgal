#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_approximation/approximate_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>                   Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;
typedef Mesh::Property_map<face_descriptor, std::size_t>      Face_proxy_pmap;

namespace VSA = CGAL::Surface_mesh_approximation;

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/bear.off");

  // read input surface triangle mesh
  Mesh mesh;
  if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, mesh) ||
     !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  // output indexed triangle mesh
  std::vector<Kernel::Point_3> anchors;
  std::vector<std::array<std::size_t, 3> > triangles; // triplets of indices

  // output face proxy index property map
  Face_proxy_pmap fpxmap =
    mesh.add_property_map<face_descriptor, std::size_t>("f:proxy_id", 0).first;

  // output planar proxies
  std::vector<Kernel::Vector_3> proxies;

  // free function interface with named parameters
  VSA::approximate_triangle_mesh(mesh,
    CGAL::parameters::min_error_drop(0.05). // seeding with minimum error drop
                      number_of_iterations(40). // set number of clustering iterations after seeding
                      subdivision_ratio(0.3). // set chord subdivision ratio threshold when meshing
                      face_proxy_map(fpxmap). // get face partition map
                      proxies(std::back_inserter(proxies)). // output proxies
                      anchors(std::back_inserter(anchors)). // output anchor points
                      triangles(std::back_inserter(triangles))); // output indexed triangles

  return EXIT_SUCCESS;
}
