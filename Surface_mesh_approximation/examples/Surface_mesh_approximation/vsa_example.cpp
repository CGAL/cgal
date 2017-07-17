#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>
#include <CGAL/vsa_mesh_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(int argc, char *argv[])
{
  if (argc < 5)
    return 0;

  // create and read Polyhedron
  Polyhedron mesh;
  std::ifstream input(argv[1]);
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // create a property-map for segment-ids
  typedef std::map<Polyhedron::Facet_const_handle, std::size_t> Facet_id_map;
  Facet_id_map internal_facet_id_map;
  for (Polyhedron::Facet_const_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    internal_facet_id_map.insert(std::pair<Polyhedron::Face_const_handle, std::size_t>(fitr, std::numeric_limits<std::size_t>::max()));
  }
  boost::associative_property_map<Facet_id_map> proxy_patch_map(internal_facet_id_map);

  const std::size_t num_proxies = std::atoi(argv[3]);
  const std::size_t num_iterations = std::atoi(argv[4]);
  std::vector<int> tris;
  std::vector<Kernel::Point_3> anchor_pos;
  std::vector<Polyhedron::Vertex_handle> anchor_vtx;
  std::vector<std::vector<std::size_t> > bdrs;
  int init = std::atoi(argv[2]);
  if (init < 0 || init > 3)
    return EXIT_FAILURE;
  CGAL::vsa_mesh_approximation(init, mesh,
    num_proxies,
    num_iterations,
    proxy_patch_map,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)),
    tris,
    anchor_pos,
    anchor_vtx,
    bdrs,
    Kernel());

  return EXIT_SUCCESS;
}

