#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_approximation/approximate_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
typedef std::map<face_descriptor, std::size_t> Face_index_map;
typedef boost::associative_property_map<Face_index_map> Face_proxy_pmap;

int main()
{
  // read input surface triangle mesh
  Polyhedron input;
  std::ifstream file("data/mask.off");
  file >> input;

  // output indexed triangle mesh
  std::vector<Kernel::Point_3> anchors;
  std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles; // triplets of indices

  // output face proxy index property map
  Face_index_map fidx_map;
  BOOST_FOREACH(face_descriptor f, faces(input))
    fidx_map[f] = 0;
  Face_proxy_pmap fpxmap(fidx_map);

  // output planar proxies
  std::vector<Kernel::Vector_3> proxies;

  // free function interface with named parameters
  CGAL::VSA::approximate_mesh(input,
    CGAL::parameters::min_error_drop(0.05). // seeding with minimum error drop
    number_of_iterations(40). // set number of clustering iterations after seeding
    subdivision_ratio(0.3). // set chord subdivision ratio threshold when meshing
    face_proxy_map(fpxmap). // get face partition map
    proxies(std::back_inserter(proxies)). // output proxies
    anchors(std::back_inserter(anchors)). // output anchor points
    triangles(std::back_inserter(triangles))); // output indexed triangles

  return EXIT_SUCCESS;
}
