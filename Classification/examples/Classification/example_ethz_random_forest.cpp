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
typedef Point_set::Property_map<unsigned char> UCmap;

namespace Classification = CGAL::Classification;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef Classification::Point_set_feature_generator<Kernel, Point_set, Pmap>    Feature_generator;


int main (int argc, char** argv)
{
  std::string filename = "data/b9_training.ply";

  if (argc > 1)
    filename = argv[1];

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

  Feature_set features;

  std::cerr << "Generating features" << std::endl;
  CGAL::Real_timer t;
  t.start();
  Feature_generator generator (pts, pts.point_map(),
                               5);  // using 5 scales

  features.begin_parallel_additions();
  generator.generate_point_based_features (features);
  features.end_parallel_additions();

  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;

  // Add labels
  Label_set labels;
  Label_handle ground = labels.add ("ground");
  Label_handle vegetation = labels.add ("vegetation");
  Label_handle roof = labels.add ("roof");

  // Check if ground truth is valid for this label set
  if (!labels.is_valid_ground_truth (pts.range(label_map), true))
    return EXIT_FAILURE;

  std::vector<int> label_indices(pts.size(), -1);

  std::cerr << "Using ETHZ Random Forest Classifier" << std::endl;
  Classification::ETHZ::Random_forest_classifier classifier (labels, features);

  std::cerr << "Training" << std::endl;
  t.reset();
  t.start();
  classifier.train (pts.range(label_map));
  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;

  t.reset();
  t.start();
  Classification::classify_with_graphcut<CGAL::Parallel_if_available_tag>
    (pts, pts.point_map(), labels, classifier,
     generator.neighborhood().k_neighbor_query(12),
     0.2f, 1, label_indices);
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

  // Color point set according to class
  UCmap red = pts.add_property_map<unsigned char>("red", 0).first;
  UCmap green = pts.add_property_map<unsigned char>("green", 0).first;
  UCmap blue = pts.add_property_map<unsigned char>("blue", 0).first;

  for (std::size_t i = 0; i < label_indices.size(); ++ i)
  {
    label_map[i] = label_indices[i]; // update label map with computed classification

    Label_handle label = labels[label_indices[i]];
    const CGAL::Color& color = label->color();
    red[i] = color.red();
    green[i] = color.green();
    blue[i] = color.blue();
  }

  // Save configuration for later use
  std::ofstream fconfig ("ethz_random_forest.bin", std::ios_base::binary);
  classifier.save_configuration(fconfig);

  // Write result
  std::ofstream f ("classification.ply");
  f.precision(18);
  f << pts;

  std::cerr << "All done" << std::endl;

  return EXIT_SUCCESS;
}
