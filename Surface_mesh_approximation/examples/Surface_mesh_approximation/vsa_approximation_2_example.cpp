#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Approximation_l21_traits<Polyhedron>::Proxy Plane_proxy;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
typedef std::map<face_descriptor, std::size_t> Facet_index_map;
typedef boost::associative_property_map<Facet_index_map> Facet_proxy_pmap;

int main()
{
  // read input surface triangle mesh 
  Polyhedron input;
  std::ifstream file("data/mask.off");
  if (!file || !(file >> input) || input.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // output indexed triangle mesh
  std::vector<Kernel::Point_3> points;
  std::vector<std::vector<std::size_t> > triangles; // triplets of indices

  // output facet proxy index property map
  Facet_index_map fidx_map;
  BOOST_FOREACH(face_descriptor f, faces(input))
    fidx_map[f] = 0;
  Facet_proxy_pmap fpxmap(fidx_map);

  // output planar proxies
  std::vector<Plane_proxy> proxies;

  // free function interface with named parameters
  CGAL::mesh_approximation(input,
    CGAL::Surface_mesh_approximation::parameters::min_error_drop(0.05). // seeding with minimum error drop
    nb_of_iterations(40). // set number of clustering iterations after seeding
    mesh_chord_error(0.3). // set chord approximation error threshold when meshing
    facet_proxy_map(fpxmap). // get facet partition map
    proxies(std::back_inserter(proxies)). // get proxies
    anchor_points(std::back_inserter(points)). // anchor points
    indexed_triangles(std::back_inserter(triangles))); // indexed triangles

  return EXIT_SUCCESS;
}
