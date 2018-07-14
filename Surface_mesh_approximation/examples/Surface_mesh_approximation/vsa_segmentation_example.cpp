#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/approximate_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::graph_traits<Polyhedron>::face_descriptor Face_descriptor;
typedef boost::unordered_map<Face_descriptor, std::size_t> Facet_index_map;
typedef boost::associative_property_map<Facet_index_map> Facet_proxy_pmap;

int main()
{
  // reads input mesh 
  Polyhedron input;
  std::ifstream file("data/mask.off");
  file >> input;

  Facet_index_map fidx_map;
  BOOST_FOREACH(Face_descriptor f, faces(input))
    fidx_map[f] = 0;

  // facet proxy index property map
  Facet_proxy_pmap fpxmap(fidx_map);

  // free function interface with named parameters
  CGAL::VSA::approximate_mesh(input,
  CGAL::VSA::parameters::max_nb_of_proxies(200). // first stop criterion
    min_error_drop(0.05). // second stop criterion
    nb_of_iterations(30). // number of relaxation iterations after seeding
    facet_proxy_map(fpxmap)); // output facet-proxy map

  // iterates over facets and outputs segment id to console
  BOOST_FOREACH(Face_descriptor f, faces(input))
    std::cout << fpxmap[f] << std::endl;

  return EXIT_SUCCESS;
}
