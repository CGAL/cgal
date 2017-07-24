#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/property_map.h>
#include <CGAL/vsa_mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_const_handle Facet_const_handle;
typedef Polyhedron::Facet_const_iterator Facet_const_iterator;
typedef std::map<Facet_const_handle, std::size_t> Facet_id_map;

int main(int argc, char *argv[])
{
  if (argc < 5)
    return 0;

  // create and read Polyhedron
  Polyhedron mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // create a property-map for facet proxy index map
  Facet_id_map facet_proxy_map;
  for (Facet_const_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr)
    facet_proxy_map.insert(std::pair<Facet_const_handle, std::size_t>(fitr, 0));
  boost::associative_property_map<Facet_id_map> f_proxy_pmap(facet_proxy_map);

  const std::size_t num_proxies = std::atoi(argv[3]);
  const std::size_t num_iterations = std::atoi(argv[4]);
  std::vector<int> tris;
  std::vector<Kernel::Point_3> anchor_pos;
  int init = std::atoi(argv[2]);
  if (init < 0 || init > 3)
    return EXIT_FAILURE;

  CGAL::vsa_approximate_and_extract(mesh,
    f_proxy_pmap,
    tris,
    anchor_pos,
    init,
    num_proxies,
    num_iterations);

  return EXIT_SUCCESS;
}
