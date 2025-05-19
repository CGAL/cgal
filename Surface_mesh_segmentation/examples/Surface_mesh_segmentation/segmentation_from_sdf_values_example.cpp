#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

int main(int argc, char* argv[])
{
    // create and read Polyhedron
    Polyhedron mesh;
    std::ifstream input(CGAL::data_file_path("meshes/cactus.off"));
    if ( !input || !(input >> mesh) || mesh.empty()  || ( !CGAL::is_triangle_mesh(mesh)))
    {
      std::cerr << "Input is not a triangle mesh." << std::endl;
      return EXIT_FAILURE;
    }

    // create a property-map for SDF values
    typedef std::map<face_descriptor, double> Facet_double_map;
    Facet_double_map internal_sdf_map;
    boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);

  // create a property-map for SDF values
  typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
  Facet_double_map internal_sdf_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);

    // create a property-map for segment-ids
    typedef std::map<face_descriptor, std::size_t> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

  // create a property-map for segment-ids
  typedef std::map<Polyhedron::Facet_const_handle, std::size_t> Facet_int_map;
  Facet_int_map internal_segment_map;
  boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

    std::cout << "Number of segments: " << number_of_segments << std::endl;
    // print segment-ids
    for(face_descriptor f : faces(mesh)) {
        // ids are between [0, number_of_segments -1]
        std::cout << segment_property_map[f] << " ";
    }
    std::cout << std::endl;

  std::cout << "Number of segments: " << number_of_segments << std::endl;
  // print segment-ids
  for(Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin();
    facet_it != mesh.facets_end(); ++facet_it) {
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
