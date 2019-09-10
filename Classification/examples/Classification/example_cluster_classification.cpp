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
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/Shape_detection_3.h>

#include <CGAL/Real_timer.h>

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

typedef Point_set::Point_map Pmap;
typedef Point_set::Vector_map Vmap;
typedef Point_set::Property_map<int> Imap;
typedef Point_set::Property_map<unsigned char> UCmap;

typedef CGAL::Shape_detection_3::Shape_detection_traits<Kernel, Point_set, Pmap, Vmap> SD_traits;
typedef CGAL::Shape_detection_3::Region_growing<SD_traits>                             Region_growing;
typedef CGAL::Shape_detection_3::Plane<SD_traits>                                      Plane;

namespace Classification = CGAL::Classification;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef Classification::Local_eigen_analysis                                    Local_eigen_analysis;
typedef Classification::Point_set_feature_generator<Kernel, Point_set, Pmap>    Feature_generator;

typedef Classification::Cluster<Point_set, Pmap>                                Cluster;


int main (int argc, char** argv)
{
  std::string filename = "data/b9.ply";
  std::string filename_config = "data/b9_clusters_config.gz";
  
  if (argc > 1)
    filename = argv[1];
  if (argc > 2)
    filename_config = argv[2];

  std::ifstream in (filename.c_str(), std::ios::binary);
  Point_set pts;

  std::cerr << "Reading input" << std::endl;
  in >> pts;

  std::cerr << "Estimating normals" << std::endl;
  CGAL::Real_timer t;
  t.start();
  pts.add_normal_map();
  CGAL::jet_estimate_normals<Concurrency_tag> (pts, 12);
  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;
  t.reset();

  Feature_set pointwise_features;
  
  std::cerr << "Generating pointwise features" << std::endl;
  t.start();
  Feature_generator generator (pts, pts.point_map(),
                               5);  // using 5 scales
  
#ifdef CGAL_LINKED_WITH_TBB
  pointwise_features.begin_parallel_additions();
#endif
  
  generator.generate_point_based_features (pointwise_features);
  generator.generate_normal_based_features (pointwise_features, pts.normal_map());

#ifdef CGAL_LINKED_WITH_TBB
  pointwise_features.end_parallel_additions();
#endif
  
  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;

  ///////////////////////////////////////////////////////////////////
  //! [Cluster]
  
  std::cerr << "Detecting planes" << std::endl;
  t.start();
  Region_growing::Parameters parameters;
  parameters.min_points = 1;
  parameters.epsilon = 1.0;
  parameters.cluster_epsilon = 1.0;
  parameters.normal_threshold = 0.9;

  Region_growing region_growing;
  region_growing.set_input (pts, pts.point_map(), pts.normal_map());
  region_growing.add_shape_factory<Plane>();
  region_growing.detect (parameters);
  t.stop();
  std::cerr << region_growing.shapes().end() - region_growing.shapes().begin() << " planes detected in "
            << t.time() << " second(s)" << std::endl;
  t.reset();

  std::cerr << "Creating clusters" << std::endl;
  t.start();
  std::vector<Cluster> clusters;
  Classification::create_clusters_from_indices (pts, pts.point_map(),
                                                CGAL::Shape_detection_3::Point_to_shape_index_map<SD_traits>
                                                (pts, region_growing.planes()),
                                                clusters);
  t.stop();
  std::cerr << clusters.size() << " clusters created in "
            << t.time() << " second(s)" << std::endl;
  t.reset();

  //! [Cluster]
  ///////////////////////////////////////////////////////////////////
  
  std::cerr << "Computing cluster features" << std::endl;
  
  ///////////////////////////////////////////////////////////////////
  //! [Eigen]
  
  Local_eigen_analysis eigen = Local_eigen_analysis::create_from_point_clusters (clusters);

  //! [Eigen]
  ///////////////////////////////////////////////////////////////////

  t.start();
  
  ///////////////////////////////////////////////////////////////////
  //! [Features]
  
  Feature_set features;
  
#ifdef CGAL_LINKED_WITH_TBB
  features.begin_parallel_additions();
#endif

  // First compute means of features
  for (std::size_t i = 0; i < pointwise_features.size(); ++ i)
    features.add<Classification::Feature::Cluster_mean_of_feature> (clusters, pointwise_features[i]);

#ifdef CGAL_LINKED_WITH_TBB
  features.end_parallel_additions();
  features.begin_parallel_additions();
#endif

  // Then compute variances of features (and remaining cluster features)
  for (std::size_t i = 0; i < pointwise_features.size(); ++ i)
    features.add<Classification::Feature::Cluster_variance_of_feature> (clusters,
                                                                        pointwise_features[i], // i^th feature
                                                                        features[i]);          // mean of i^th feature

  features.add<Classification::Feature::Cluster_size> (clusters);
  features.add<Classification::Feature::Cluster_vertical_extent> (clusters);
  
  for (std::size_t i = 0; i < 3; ++ i)
    features.add<Classification::Feature::Eigenvalue> (clusters, eigen, (unsigned int)(i));
  
#ifdef CGAL_LINKED_WITH_TBB
  features.end_parallel_additions();
#endif
  
  //! [Features]
  ///////////////////////////////////////////////////////////////////
  
  t.stop();
  
  // Add types
  Label_set labels;
  Label_handle ground = labels.add ("ground");
  Label_handle vegetation = labels.add ("vegetation");
  Label_handle roof = labels.add ("roof");

  std::vector<int> label_indices(clusters.size(), -1);
  
  std::cerr << "Using ETHZ Random Forest Classifier" << std::endl;
  Classification::ETHZ_random_forest_classifier classifier (labels, features);
  
  std::cerr << "Loading configuration" << std::endl;
  std::ifstream in_config (filename_config, std::ios_base::in | std::ios_base::binary);
  classifier.load_configuration (in_config);

  std::cerr << "Classifying" << std::endl;
  t.reset();
  t.start();
  Classification::classify<Concurrency_tag> (clusters, labels, classifier, label_indices);
  t.stop();
  
  std::cerr << "Classification done in " << t.time() << " second(s)" << std::endl;
  
  return EXIT_SUCCESS;
}
