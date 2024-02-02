#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

#define TEST_CLASSIFIER false

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <bitset>

#define CGAL_CLASSIFICATION_VERBOSE
#define CGAL_COVERAGE_FEATURE_DISPLAY_PROGRESS true

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Random.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Point_map Point_map;

typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

namespace Classification = CGAL::Classification;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef Classification::Sum_of_weighted_features_classifier                     Classifier;

typedef Classification::Point_set_feature_generator<Kernel, Point_set, Point_map> Feature_generator;

typedef Point_set::Vector_map Vector_map;
typedef Point_set::Property_map<std::size_t> Size_t_map;
typedef Point_set::Property_map<CGAL::IO::Color> Color_map;

typedef Point_set::Point_map Pmap;
typedef Point_set::Property_map<int> Imap;
typedef Point_set::Property_map<float> Fmap;
typedef Point_set::Property_map<unsigned char> UCmap;


#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
namespace AW3 = CGAL::Alpha_wraps_3;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_Alpha;
typedef Kernel_Alpha::Point_3 Point_Alpha;
using Point_container = std::vector<Point_Alpha>;
using Mesh = CGAL::Surface_mesh<Point_Alpha>;


int main (int argc, char** argv) //(int, char**)
{
  //const std::string filename_train = (argc > 1) ? argv[1] : "C:/Users/rdyke/Documents/Datasets/Semantic3D/bildstein_station1_xyz_intensity_rgb - crop - small int.ply";
  const std::string filename_train = (argc > 1) ? argv[1] : "C:/dev/Point-cloud-descriptors-for-classification/Data/data/points_3/b9_training.ply";
  //const std::string filename_train = (argc > 1) ? argv[1] : "C:/Users/rdyke/Documents/Datasets/Semantic3D/b9_training_shifted.ply";
  //const std::string filename_train = "C:/Users/rdyke/Documents/Datasets/Semantic3D/bildstein_station1_xyz_intensity_rgb - smaller.ply";
  const std::string filename_test = (argc > 2) ? argv[2] : "C:/Users/rdyke/Documents/Datasets/Semantic3D/bildstein_station3_xyz_intensity_rgb - crop - small int.ply";
  const std::string filename_results = (argc > 3) ? argv[3] : "results.txt";
  //const int options = (argc > 4) ? atoi(argv[4]) : 0b0000011;
  //const int options = (argc > 4) ? atoi(argv[4]) : 0b0100000;
  const int options = (argc > 4) ? atoi(argv[4]) : 0b0100000;
  // options:
  //  point-based features      0b0000001
  //  color-based features      0b0000010
  //  betti numbers             0b0000100
  //  compactness               0b0001000
  //  coverage                  0b0010000
  //  fractal dimensionality    0b0100000
  //  skeletonisation           0b1000000
  std::bitset<7> boptions(options);

  std::cout << "options:" << std::endl;
  std::cout << "\tpoint-based features " << ((boptions[0]) ? "enabled" : "disabled") << std::endl;
  std::cout << "\tcolor-based features " << ((boptions[1]) ? "enabled" : "disabled") << std::endl;
  std::cout << "\tbetti numbers " << ((boptions[2]) ? "enabled" : "disabled") << std::endl;
  std::cout << "\tcompactness " << ((boptions[3]) ? "enabled" : "disabled") << std::endl;
  std::cout << "\tcoverage " << ((boptions[4]) ? "enabled" : "disabled") << std::endl;
  std::cout << "\tfractal dimensionality " << ((boptions[5]) ? "enabled" : "disabled") << std::endl;
  std::cout << "\tskeletonisation " << ((boptions[6]) ? "enabled" : "disabled") << std::endl;
  if (options == 0) {
      std::cerr << "No features selected" << std::endl;
      return EXIT_FAILURE;
  }

  std::ofstream f_results(filename_results);

  f_results << filename_train << std::endl << filename_test << std::endl;

  std::ifstream in(filename_train.c_str(), std::ios::binary);
  if (!in) {
    std::cerr << "Failed to open file." << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Timer t;
  t.reset(); t.start();

  Point_set pts;

  in >> pts;

  assert(pts.number_of_points() > 0);
  std::cout << "size = " << pts.number_of_points() << std::endl;

  Imap label_map;
  bool lm_found = false;
  std::tie(label_map, lm_found) = pts.property_map<int>("label");
  //assert(lm_found);
  if (!lm_found) {
    std::cerr << "Error: Ground truth labels not found." << std::endl;
    return EXIT_FAILURE;
  }

  Color_map color_map;
  bool cm_found = false;
  boost::tie(color_map, cm_found) = pts.add_property_map<CGAL::IO::Color>("color");
  //assert(color_map);


  Mesh wrap;
  if (boptions[2] || boptions[3] || boptions[4] || boptions[6]) {
      // wrap point cloud
      const float radius_size = 0.1f; // specify radius of neighborhoods (default: 10cm==0.1f)
      const float voxel_size = radius_size / 3.f; // re-scale for CGAL's feature generator
      const double relative_alpha = 0.125f; //0.2f
      const double relative_offset = 2.f; //2.f
      double alpha = radius_size / relative_alpha;
      double offset = radius_size / relative_offset;

      // convert to a kernel that is more stable for Alpha Wrap
      CGAL::Cartesian_converter<Kernel, Kernel_Alpha> to_mesh;
      Point_container points;
      points.reserve(pts.number_of_points());
      for (const auto& point : pts.points()) {
          points.push_back(to_mesh(point));
      }

      CGAL::alpha_wrap_3(points, alpha, offset, wrap);
      std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces, " << std::endl;

      CGAL::IO::write_polygon_mesh("wrap.ply", wrap, CGAL::parameters::stream_precision(17));
      std::cout << "Wrap saved" << std::endl;
  }


  Feature_set features;

  Feature_generator generator(pts, pts.point_map(), 5);  // using 5 scales

#ifdef CGAL_LINKED_WITH_TBB
  features.begin_parallel_additions();
#endif

  if (boptions[0])
    generator.generate_point_based_features(features);
  if (boptions[1])
    generator.generate_color_based_features(features, color_map);

#ifdef CGAL_LINKED_WITH_TBB
  features.end_parallel_additions();
#endif

  if (boptions[2]) {
      std::cout << "Computing betti feature..." << std::endl;
      using MyBettiNumbers = CGAL::Classification::Feature::Betti_numbers<Kernel, Point_set, Point_set::Point_map, Mesh>;
      auto bbox = CGAL::bounding_box(CGAL::make_transform_iterator_from_property_map(pts.begin(), pts.point_map()),
          CGAL::make_transform_iterator_from_property_map(pts.end(), pts.point_map()));
      for (std::size_t i = 0; i < generator.number_of_scales(); ++i) {
          std::cout << "  scale: " << i << std::endl;
          features.add_multidimensional_feature_with_scale_id<MyBettiNumbers>(i, 2, pts, pts.point_map(), generator.grid(i), generator.radius_neighbors(i), wrap);
      }
  }
  if (boptions[3]) {
      std::cout << "Computing compactness feature..." << std::endl;
      using MyCompactness = CGAL::Classification::Feature::Compactness<Kernel, Point_set, Point_set::Point_map, Mesh>;
      auto bbox = CGAL::bounding_box(CGAL::make_transform_iterator_from_property_map(pts.begin(), pts.point_map()),
          CGAL::make_transform_iterator_from_property_map(pts.end(), pts.point_map()));
      using Planimetric_grid = Classification::Planimetric_grid<Kernel, Point_set, Pmap>;
      for (std::size_t i = 0; i < generator.number_of_scales(); ++i) {
          std::cout << "  scale: " << i << std::endl;
          features.add_with_scale_id<MyCompactness>(i, pts, pts.point_map(), generator.grid(i), generator.radius_neighbors(i), wrap);
      }
  }
  if (boptions[4]) {
      std::cout << "Computing coverage feature..." << std::endl;
      using MyCoverage = CGAL::Classification::Feature::Coverage<Kernel, Point_set, Pmap, Mesh>;
      for (std::size_t i = 0; i < generator.number_of_scales(); ++i) {
          std::cout << "  scale: " << i << std::endl;
          features.add_with_scale_id<MyCoverage>(i, pts, pts.point_map(), wrap, generator.radius_neighbors(i));
      }
  }
  if (boptions[5]) {
      std::cout << "Computing fractal feature..." << std::endl;
      using MyFractalDimensionality = CGAL::Classification::Feature::Fractal_dimensionality<Kernel, Point_set, Point_set::Point_map>;
      for (std::size_t i = 0; i < generator.number_of_scales(); ++i) {
          features.add_with_scale_id<MyFractalDimensionality>(i, pts, pts.point_map(), generator.radius_neighbors(i));
      }
  }
  if (boptions[6]) {
      std::cout << "Computing skeletonize feature..." << std::endl;
      using MySkeleton = CGAL::Classification::Feature::Skeleton<Kernel, Point_set, Pmap, Mesh>;
      features.add_multidimensional_feature<MySkeleton>(2, pts, pts.point_map(), wrap);
  }
  std::cout << "Features computed" << std::endl;
  
  assert (generator.number_of_scales() == 5);
  //assert (features.size() == 59);
  std::cout << "features.size() = " << features.size() << std::endl;



  // Export features
  std::cout << "pts.number_of_points() = " << pts.number_of_points() << ", pts.size() = " << pts.size() << std::endl;
  for (auto feature : features) {
      std::cout << "feature name: " << feature->name() << std::endl;
      Fmap prop = pts.add_property_map<float>("scalar_" + feature->name(), 0).first;
      for (size_t i = 0; i < pts.number_of_points(); ++i) {
          prop[i] = feature->value(i);
      }
  }

  pts.remove_property_map(pts.property_map<unsigned char>("red").first);
  pts.remove_property_map(pts.property_map<unsigned char>("green").first);
  pts.remove_property_map(pts.property_map<unsigned char>("blue").first);


  std::ofstream f_train2("features_train.ply");
  f_train2.precision(12);
  f_train2 << pts;

  t.stop();
  std::cout << "Export stage complete (" << t.time() << "s)" << std::endl;


  return EXIT_SUCCESS;

  Label_set labels;
  Label_handle undefined = labels.add("undefined", CGAL::IO::Color(0, 0, 0), 0);
  Label_handle manmade_terrain = labels.add("man-made terrain", CGAL::IO::Color(0, 0, 255), 1);
  Label_handle natural_terrain = labels.add("natural terrain", CGAL::IO::Color(0, 125, 255), 2);
  Label_handle high_vegetation = labels.add("high vegetation", CGAL::IO::Color(125, 125, 0), 3);
  Label_handle low_vegetation = labels.add("low vegetation", CGAL::IO::Color(0, 255, 0), 4);
  Label_handle buildings = labels.add("buildings", CGAL::IO::Color(255, 165, 0), 5);
  Label_handle hard_scape = labels.add("hard scape", CGAL::IO::Color(125, 165, 125), 6);
  Label_handle scanning_artefacts = labels.add("scanning artefacts", CGAL::IO::Color(255, 0, 0), 7);
  Label_handle cars = labels.add("cars", CGAL::IO::Color(165, 165, 165), 8);

  if (!labels.is_valid_ground_truth(pts.range(label_map), true)) {
      std::cerr << "Error: Ground truth labels are invalid." << std::endl;
      return EXIT_FAILURE;
  }

  //assert (labels.size() == 20);
  std::cout << "labels.size() = " << labels.size() << std::endl;

  std::vector<int> label_indices(pts.size(), -1);

  //Classifier classifier (labels, features);
  Classification::ETHZ::Random_forest_classifier classifier(labels, features);

  classifier.train<CGAL::Parallel_if_available_tag> (pts.range(label_map));

  std::cout << "Feature usage:" << std::endl;
  std::vector<std::size_t> count;
  classifier.get_feature_usage(count);
  int total = 0;
  for (auto c : count)
    total += c;
  for (int i = 0; i < features.size(); ++i) {
    std::cout << "\t" << features[i]->name() << ": " << count[i] << "/" << total << " (" << 100.f * count[i] / total << "%)" << std::endl;
  }

  Classification::classify<CGAL::Parallel_if_available_tag>
    (pts, labels, classifier, label_indices);
  /*
  Classification::classify_with_local_smoothing<CGAL::Parallel_if_available_tag>
    (pts, pts.point_map(), labels, classifier,
     generator.neighborhood().sphere_neighbor_query(0.01f),
     label_indices);

  Classification::classify_with_graphcut<CGAL::Parallel_if_available_tag>
    (pts, pts.point_map(), labels, classifier,
     generator.neighborhood().k_neighbor_query(12),
     0.2f, 10, label_indices);
  */

  // Save configuration for later use
  std::ofstream fconfig("ethz_random_forest.bin", std::ios_base::binary);
  classifier.save_configuration(fconfig);
  fconfig.close();

  Classification::Evaluation evaluation (labels, pts.range(label_map), label_indices);
  
  for (Label_handle l : labels) {
    std::cout << " * " << l->name() << ": "
              << evaluation.precision(l) << " ; "
              << evaluation.recall(l) << " ; "
              << evaluation.f1_score(l) << " ; "
              << evaluation.intersection_over_union(l) << std::endl;
  }

  std::cout << evaluation << std::endl;
  std::cout << "(train) Accuracy = " << evaluation.accuracy() << std::endl
            << "(train) Mean F1 score = " << evaluation.mean_f1_score() << std::endl
            << "(train) Mean IoU = " << evaluation.mean_intersection_over_union() << std::endl;

  f_results << "(train) accuracy," << evaluation.accuracy() << std::endl
            << "(train) Mean F1 score," << evaluation.mean_f1_score() << std::endl
            << "(train) Mean IoU," << evaluation.mean_intersection_over_union() << std::endl;

  Imap label_pred_map = pts.add_property_map<int>("labels_pred", 0).first;
  Imap label_valid_map = pts.add_property_map<int>("label_valid", 0).first;
  for (std::size_t i = 0; i < label_indices.size(); ++i) {
    label_pred_map[i] = label_indices[i];
    label_valid_map[i] = (label_indices[i] == label_map[i]) ? 1 : (label_map[i] == -1) ? -1 : 0;
  }


  pts.remove_property_map(pts.property_map<unsigned char>("red").first);
  pts.remove_property_map(pts.property_map<unsigned char>("green").first);
  pts.remove_property_map(pts.property_map<unsigned char>("blue").first);

  
  std::ofstream f_train("classification_train.ply");
  f_train.precision(12);
  f_train << pts;
  
  t.stop();
  std::cout << "Training stage complete (" << t.time() << "s)" << std::endl;


#if TEST_CLASSIFIER
  pts.clear();

  //in.open(filename_test.c_str(), std::ios::binary);
  std::ifstream in_test(filename_test.c_str(), std::ios::binary);
  if (!in_test) {
    std::cerr << "Failed to open file." << std::endl;
    return EXIT_FAILURE;
  }
  in_test >> pts;

  std::cout << "size = " << pts.number_of_points() << std::endl;

  std::tie(label_map, lm_found) = pts.property_map<int>("label");
  //assert(lm_found);

  boost::tie(color_map, cm_found) = pts.add_property_map<CGAL::IO::Color>("color");
  //assert(color_map);

  Feature_set features_test;

  Feature_generator generator_test(pts, pts.point_map(), 5);  // using 5 scales

#ifdef CGAL_LINKED_WITH_TBB
  features_test.begin_parallel_additions();
#endif

  generator_test.generate_point_based_features(features_test);
  generator_test.generate_color_based_features(features_test, color_map);

#ifdef CGAL_LINKED_WITH_TBB
  features_test.end_parallel_additions();
#endif

  std::vector<int> label_test_indices(pts.size(), -1);

  Classification::ETHZ::Random_forest_classifier classifier_test(labels, features_test);

  std::ifstream in_fconfig("ethz_random_forest.bin", std::ios_base::binary);
  classifier_test.load_configuration(in_fconfig);
  in_fconfig.close();

  Classification::classify<CGAL::Parallel_if_available_tag>
      (pts, labels, classifier_test, label_test_indices);

  Classification::Evaluation evaluation_test (labels, pts.range(label_map), label_test_indices);

  for (Label_handle l : labels) {
    std::cout << " * " << l->name() << ": "
              << evaluation_test.precision(l) << " ; "
              << evaluation_test.recall(l) << " ; "
              << evaluation_test.f1_score(l) << " ; "
              << evaluation_test.intersection_over_union(l) << std::endl;
  }

  std::cout << evaluation_test << std::endl;
  std::cout << "(test) Accuracy = " << evaluation_test.accuracy() << std::endl
            << "(test) Mean F1 score = " << evaluation_test.mean_f1_score() << std::endl
            << "(test) Mean IoU = " << evaluation_test.mean_intersection_over_union() << std::endl;

  f_results << "(test) accuracy," << evaluation_test.accuracy() << std::endl
            << "(test) Mean F1 score," << evaluation_test.mean_f1_score() << std::endl
            << "(test) Mean IoU," << evaluation_test.mean_intersection_over_union() << std::endl;

  
  for (int i = 0; i < features.size(); ++i) { // feature usage from training
      f_results << "(feature usage) " << features[i]->name() << "," << count[i] << "," << total << std::endl;
  }
  f_results.close();

  label_pred_map = pts.add_property_map<int>("labels_pred", 0).first;
  label_valid_map = pts.add_property_map<int>("label_valid", 0).first;
  for (std::size_t i = 0; i < label_test_indices.size(); ++i) {
    label_pred_map[i] = label_test_indices[i];
    label_valid_map[i] = (label_test_indices[i] == label_map[i]) ? 1 : (label_map[i] == -1) ? -1 : 0;
  }

  pts.remove_property_map(pts.property_map<unsigned char>("red").first);
  pts.remove_property_map(pts.property_map<unsigned char>("green").first);
  pts.remove_property_map(pts.property_map<unsigned char>("blue").first);

  std::ofstream f_test("classification_test.ply");
  f_test.precision(12);
  f_test << pts;
#endif
  return EXIT_SUCCESS;
}