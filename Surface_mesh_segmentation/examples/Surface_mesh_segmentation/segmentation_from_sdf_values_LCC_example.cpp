#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/mesh_segmentation.h>

#include <CGAL/property_map.h>
#include <CGAL/Unique_hash_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
    <2, 3, MyTraits>::type LCC;
typedef LCC::Attribute_const_handle<2>::type Facet_const_handle;

int main()
{
  // create and read LCC
  LCC mesh;
  CGAL::read_off("data/cactus.off", mesh);

  // create a property-map for SDF values
  typedef CGAL::Unique_hash_map<Facet_const_handle, double> Facet_double_map;
  Facet_double_map internal_sdf_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);

  // compute SDF values using default parameters for number of rays, and cone angle
  CGAL::sdf_values(mesh, sdf_property_map);

  // create a property-map for segment-ids
  typedef CGAL::Unique_hash_map<Facet_const_handle, std::size_t> Facet_int_map;
  Facet_int_map internal_segment_map;
  boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

  // segment the mesh using default parameters for number of levels, and smoothing lambda
  // Any other scalar values can be used instead of using SDF values computed using the CGAL function
  std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map,
                                    segment_property_map);

  std::cout << "Number of segments: " << number_of_segments << std::endl;
  // print segment-ids
  for(LCC::Attribute_range<2>::type::const_iterator facet_it=mesh.attributes<2>().begin();
    facet_it!=mesh.attributes<2>().end(); ++facet_it)
  {
   // ids are between [0, number_of_segments -1]
   std::cout << segment_property_map[facet_it] << " ";
  }
  std::cout << std::endl;

  const std::size_t number_of_clusters = 4;       // use 4 clusters in soft clustering
  const double smoothing_lambda = 0.3;  // importance of surface features, suggested to be in-between [0,1]

  // Note that we can use the same SDF values (sdf_property_map) over and over again for segmentation.
  // This feature is relevant for segmenting the mesh several times with different parameters.
  CGAL::segmentation_from_sdf_values(
   mesh, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda);

  return EXIT_SUCCESS;
}
