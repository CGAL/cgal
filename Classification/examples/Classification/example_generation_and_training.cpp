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

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

typedef Point_set::Point_map Pmap;
typedef Point_set::Property_map<int> Imap;

namespace Classification = CGAL::Classification;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef Classification::Sum_of_weighted_features_classifier                     Classifier;

typedef Classification::Point_set_feature_generator<Kernel, Point_set, Pmap>    Feature_generator;


int main (int argc, char** argv)
{
  std::string filename (argc > 1 ? argv[1] : "data/b9_training.ply");
  std::ifstream in (filename.c_str(), std::ios::binary);
  Point_set pts;

  std::cerr << "Reading input" << std::endl;
  in >> pts;

  Imap label_map;
  bool lm_found = false;
  std::tie (label_map, lm_found) = pts.property_map<int> ("label");
  if (!lm_found)
  {
    std::cerr << "Error: \"label\" property not found in input file." << std::endl;
    return EXIT_FAILURE;
  }

  std::cerr << "Generating features" << std::endl;
  CGAL::Real_timer t;
  t.start();

  ///////////////////////////////////////////////////////////////////
  //! [Generator]

  Feature_set features;

  std::size_t number_of_scales = 5;
  Feature_generator generator (pts, pts.point_map(), number_of_scales);

  features.begin_parallel_additions();
  generator.generate_point_based_features (features);
  features.end_parallel_additions();

  //! [Generator]
  ///////////////////////////////////////////////////////////////////

  t.stop();
  std::cerr << features.size() << " feature(s) generated in " << t.time() << " second(s)" << std::endl;

  Label_set labels = { "ground", "vegetation", "roof" };

  Classifier classifier (labels, features);

  std::cerr << "Training" << std::endl;
  t.reset();
  t.start();
  classifier.train<CGAL::Parallel_if_available_tag> (pts.range(label_map), 800);
  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;

  t.reset();
  t.start();
  std::vector<int> label_indices(pts.size(), -1);
  Classification::classify_with_graphcut<CGAL::Parallel_if_available_tag>
    (pts, pts.point_map(), labels, classifier,
     generator.neighborhood().k_neighbor_query(12),
     0.2f, 10, label_indices);
  t.stop();
  std::cerr << "Classification with graphcut done in " << t.time() << " second(s)" << std::endl;

  std::cerr << "Precision, recall, F1 scores and IoU:" << std::endl;
  Classification::Evaluation evaluation (labels, pts.range(label_map), label_indices);

  for (Label_handle l : labels)
  {
    std::cerr << " * " << l->name() << ": "
              << evaluation.precision(l) << " ; "
              << evaluation.recall(l) << " ; "
              << evaluation.f1_score(l) << " ; "
              << evaluation.intersection_over_union(l) << std::endl;
  }

  std::cerr << "Accuracy = " << evaluation.accuracy() << std::endl
            << "Mean F1 score = " << evaluation.mean_f1_score() << std::endl
            << "Mean IoU = " << evaluation.mean_intersection_over_union() << std::endl;


  /// Save the configuration to be able to reload it later
  std::ofstream fconfig ("config.xml");
  classifier.save_configuration (fconfig);
  fconfig.close();

  std::cerr << "All done" << std::endl;

  return EXIT_SUCCESS;
}
