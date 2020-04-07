#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <boost/function_output_iterator.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/Shape_detection/Region_growing.h>
#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3                Point;
typedef Kernel::Iso_cuboid_3           Iso_cuboid_3;

typedef CGAL::Point_set_3<Point> Point_set;

typedef Point_set::Point_map                   Pmap;
typedef Point_set::Vector_map                  Vmap;
typedef Point_set::Property_map<int>           Imap;
typedef Point_set::Property_map<unsigned char> UCmap;


typedef CGAL::Shape_detection::Point_set::Sphere_neighbor_query<Kernel, Point_set, Pmap>                Neighbor_query;
typedef CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region<Kernel, Point_set, Pmap, Vmap> Region_type;
typedef CGAL::Shape_detection::Region_growing<Point_set, Neighbor_query, Region_type>                   Region_growing;

namespace Classification = CGAL::Classification;
namespace Feature = CGAL::Classification::Feature;

typedef Classification::Label_handle   Label_handle;
typedef Classification::Feature_handle Feature_handle;
typedef Classification::Label_set      Label_set;
typedef Classification::Feature_set    Feature_set;

typedef Classification::Local_eigen_analysis                                 Local_eigen_analysis;
typedef Classification::Point_set_feature_generator<Kernel, Point_set, Pmap> Feature_generator;
typedef Classification::Cluster<Point_set, Pmap>                             Cluster;

int main (int argc, char** argv)
{
  std::string filename        = "data/b9.ply";
  std::string filename_config = "data/b9_clusters_config.bin";

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
  CGAL::jet_estimate_normals<CGAL::Parallel_if_available_tag> (pts, 12);
  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;
  t.reset();

  Feature_set pointwise_features;

  std::cerr << "Generating pointwise features" << std::endl;
  t.start();
  Feature_generator generator (pts, pts.point_map(), 5); // using 5 scales

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

  std::cerr << "Detecting planes and creating clusters" << std::endl;
  t.start();

  const double search_sphere_radius  = 1.0;
  const double max_distance_to_plane = 1.0;
  const double max_accepted_angle    = 25.0;
  const std::size_t min_region_size  = 10;

  Neighbor_query neighbor_query (
    pts,
    search_sphere_radius,
    pts.point_map());
  Region_type region_type (
    pts,
    max_distance_to_plane, max_accepted_angle, min_region_size,
    pts.point_map(), pts.normal_map());
  Region_growing region_growing (
    pts, neighbor_query, region_type);

  std::vector<Cluster> clusters;
  region_growing.detect
    (boost::make_function_output_iterator
     ([&](const std::vector<std::size_t>& region) -> void {

        // Create a new cluster.
        Classification::Cluster<Point_set, Pmap> cluster (pts, pts.point_map());
        for (const std::size_t idx : region)
          cluster.insert(idx);
        clusters.push_back(cluster);
      }));

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

  // First, compute means of features.
  features.begin_parallel_additions();
  for (Feature_handle fh : pointwise_features)
    features.add<Feature::Cluster_mean_of_feature> (clusters, fh);
  features.end_parallel_additions();

  // Then, compute variances of features (and remaining cluster features).
  features.begin_parallel_additions();
  for (std::size_t i = 0; i < pointwise_features.size(); ++ i)
    features.add<Feature::Cluster_variance_of_feature> (clusters,
                                                        pointwise_features[i], // i^th feature
                                                        features[i]);          // mean of i^th feature

  features.add<Feature::Cluster_size> (clusters);
  features.add<Feature::Cluster_vertical_extent> (clusters);

  for (std::size_t i = 0; i < 3; ++ i)
    features.add<Feature::Eigenvalue> (clusters, eigen, (unsigned int)(i));

  features.end_parallel_additions();

  //! [Features]
  ///////////////////////////////////////////////////////////////////

  t.stop();

  Label_set labels = { "ground", "vegetation", "roof" };

  std::vector<int> label_indices(clusters.size(), -1);

  std::cerr << "Using ETHZ Random Forest Classifier" << std::endl;
  Classification::ETHZ::Random_forest_classifier classifier (labels, features);

  std::cerr << "Loading configuration" << std::endl;
  std::ifstream in_config (filename_config, std::ios_base::in | std::ios_base::binary);
  classifier.load_configuration (in_config);

  std::cerr << "Classifying" << std::endl;
  t.reset();
  t.start();
  Classification::classify<CGAL::Parallel_if_available_tag> (clusters, labels, classifier, label_indices);
  t.stop();

  std::cerr << "Classification done in " << t.time() << " second(s)" << std::endl;
  return EXIT_SUCCESS;
}
