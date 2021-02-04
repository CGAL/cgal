#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/mesh_segmentation.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point_3;

typedef CGAL::Surface_mesh<Point_3>                           Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;

int main(int argc, char** argv )
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cactus.off");

  Mesh mesh;
  if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, mesh) ||
     !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }

  typedef Mesh::Property_map<face_descriptor,double> Facet_double_map;
  Facet_double_map sdf_property_map;

  sdf_property_map = mesh.add_property_map<face_descriptor,double>("f:sdf").first;

  // compute SDF values
  // We can't use default parameters for number of rays, and cone angle
  // and the postprocessing
  CGAL::sdf_values(mesh, sdf_property_map, 2.0 / 3.0 * CGAL_PI, 25, true);

  // create a property-map for segment-ids
  typedef Mesh::Property_map<face_descriptor, std::size_t> Facet_int_map;
  Facet_int_map segment_property_map = mesh.add_property_map<face_descriptor,std::size_t>("f:sid").first;;

  // segment the mesh using default parameters for number of levels, and smoothing lambda
  // Any other scalar values can be used instead of using SDF values computed using the CGAL function
  std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);

  std::cout << "Number of segments: " << number_of_segments << std::endl;
  // print segment-ids

  for(face_descriptor fd : faces(mesh)){
      // ids are between [0, number_of_segments -1]
      std::cout << segment_property_map[fd] << " ";
  }
  std::cout << std::endl;

  const std::size_t number_of_clusters = 4;       // use 4 clusters in soft clustering
  const double smoothing_lambda = 0.3;  // importance of surface features, suggested to be in-between [0,1]

  // Note that we can use the same SDF values (sdf_property_map) over and over again for segmentation.
  // This feature is relevant for segmenting the mesh several times with different parameters.
  CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda);

  return EXIT_SUCCESS;
}
