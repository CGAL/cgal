#define CGAL_BGL_TESTSUITE
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#include <CGAL/mesh_segmentation.h>

#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Mesh;

typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

int main(int argc, char** argv )
{
  Mesh mesh;
  if (argc==2)
    OpenMesh::IO::read_mesh(mesh, argv[1]);
  else
    OpenMesh::IO::read_mesh(mesh, "data/cactus.off");

  if (!CGAL::is_triangle_mesh(mesh)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "#F : " << num_faces(mesh) << std::endl;
  std::cout << "#H : " << num_halfedges(mesh) << std::endl;
  std::cout << "#V : " << num_vertices(mesh) << std::endl;

  // create a property-map for SDF values
  typedef std::map<face_descriptor, double> Facet_double_map;
  Facet_double_map internal_sdf_map;
  boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);

  // compute SDF values
  CGAL::sdf_values(mesh, sdf_property_map);

  // create a property-map for segment-ids
  typedef std::map<face_descriptor, std::size_t> Facet_int_map;
  Facet_int_map internal_segment_map;
  boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

  // segment the mesh using default parameters for number of levels, and smoothing lambda
  // Any other scalar values can be used instead of using SDF values computed using the CGAL function
  std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);

  std::cout << "Number of segments: " << number_of_segments << std::endl;
  // print segment-ids

  for(face_descriptor f : faces(mesh)) {
      // ids are between [0, number_of_segments -1]
      std::cout << segment_property_map[f] << " ";
  }
  std::cout << std::endl;

  const std::size_t number_of_clusters = 4;       // use 4 clusters in soft clustering
  const double smoothing_lambda = 0.3;  // importance of surface features, suggested to be in-between [0,1]

  // Note that we can use the same SDF values (sdf_property_map) over and over again for segmentation.
  // This feature is relevant for segmenting the mesh several times with different parameters.
  CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda);
  return EXIT_SUCCESS;
}
