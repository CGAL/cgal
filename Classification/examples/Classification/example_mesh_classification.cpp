#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;

namespace Classification = CGAL::Classification;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef Classification::Face_descriptor_to_center_of_mass_map<Mesh>             Face_point_map;
typedef Classification::Face_descriptor_to_face_descriptor_with_bbox_map<Mesh>  Face_with_bbox_map;
typedef Classification::Mesh_feature_generator<Kernel, Mesh, Face_point_map>    Feature_generator;

int main (int argc, char** argv)
{
  std::string filename = "data/b9_mesh.off";
  std::string filename_config = "data/b9_mesh_config.bin";

  if (argc > 1)
    filename = argv[1];
  if (argc > 2)
    filename_config = argv[2];

  std::ifstream in (filename.c_str());
  Mesh mesh;

  std::cerr << "Reading input" << std::endl;
  in >> mesh;

  std::cerr << "Generating features" << std::endl;
  CGAL::Real_timer t;
  t.start();

  ///////////////////////////////////////////////////////////////////
  //! [Generator]

  Feature_set features;

  Face_point_map face_point_map (&mesh); // Associates each face to its center of mass

  std::size_t number_of_scales = 5;
  Feature_generator generator (mesh, face_point_map, number_of_scales);

  features.begin_parallel_additions();
  generator.generate_point_based_features (features); // Features that consider the mesh as a point set
  generator.generate_face_based_features (features);  // Features computed directly on mesh faces
  features.end_parallel_additions();

  //! [Generator]
  ///////////////////////////////////////////////////////////////////

  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;

  Label_set labels = { "ground", "vegetation", "roof" };

  std::vector<int> label_indices(mesh.number_of_faces(), -1);

  std::cerr << "Using ETHZ Random Forest Classifier" << std::endl;
  Classification::ETHZ::Random_forest_classifier classifier (labels, features);

  std::cerr << "Loading configuration" << std::endl;
  std::ifstream in_config (filename_config, std::ios_base::in | std::ios_base::binary);
  classifier.load_configuration (in_config);

  std::cerr << "Classifying with graphcut" << std::endl;
  t.reset();
  t.start();
  Classification::classify_with_graphcut<CGAL::Parallel_if_available_tag>
    (mesh.faces(), Face_with_bbox_map(&mesh), labels, classifier,
     generator.neighborhood().n_ring_neighbor_query(2),
     0.2f, 1, label_indices);
  t.stop();

  std::cerr << "Classification with graphcut done in " << t.time() << " second(s)" << std::endl;

  return EXIT_SUCCESS;
}
