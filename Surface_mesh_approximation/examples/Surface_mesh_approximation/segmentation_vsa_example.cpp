#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>
#include <CGAL/internal/Surface_mesh_approximation/VSA_segmentation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(int argc, char *argv[])
{
  if (argc < 4)
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
  Facet_id_map internal_segment_map;
  for (Polyhedron::Facet_const_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    internal_segment_map.insert(std::pair<Polyhedron::Face_const_handle, std::size_t>(fitr, std::numeric_limits<std::size_t>::max()));
  }
  boost::associative_property_map<Facet_id_map> segment_property_map(internal_segment_map);

  typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type PointPropertyMap;
  PointPropertyMap ppmap = get(boost::vertex_point, const_cast<Polyhedron &>(mesh));

  CGAL::internal::VSA_segmentation<Polyhedron, Kernel, PointPropertyMap> vsa_seg(mesh, ppmap, Kernel());

  /*for (Polyhedron::Facet_const_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    std::cout << segment_property_map[fitr] << '\n';
  }*/

  vsa_seg.partition(std::atoi(argv[2]), std::atoi(argv[3]), segment_property_map);

  int count = 0;
  for (Polyhedron::Facet_const_iterator fitr = mesh.facets_begin(); fitr != mesh.facets_end(); ++fitr) {
    //std::cout << segment_property_map[fitr] << '\n';
    ++count;
  }
  std::cout << count << std::endl;
}

