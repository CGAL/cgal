#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//#define CGAL_CLASSIFICATION_VERBOSE

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Classification/Random_forest_predicate.h>
#include <CGAL/IO/read_ply_points.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point> Point_range;
typedef CGAL::Identity_property_map<Point> Pmap;

namespace Classif = CGAL::Classification;

typedef Classif::Label_handle                                            Label_handle;
typedef Classif::Feature_handle                                          Feature_handle;
typedef Classif::Label_set                                               Label_set;
typedef Classif::Feature_set                                             Feature_set;

typedef Classif::Random_forest_predicate Classification_predicate;

typedef Classif::Point_set_feature_generator<Kernel, Point_range, Pmap> Feature_generator;

/*
 This interpreter is used to read a PLY input that contains training
 attributes (with the PLY "label" property).
*/
class My_ply_interpreter
{
  std::vector<Point>& points;
  std::vector<std::size_t>& labels;
    
public:
  My_ply_interpreter (std::vector<Point>& points,
                      std::vector<std::size_t>& labels)
    : points (points), labels (labels)
  { }

  // Init and test if input file contains the right properties
  bool is_applicable (CGAL::Ply_reader& reader)
  {
    return reader.does_tag_exist<double> ("x")
      && reader.does_tag_exist<double> ("y")
      && reader.does_tag_exist<double> ("z")
      && reader.does_tag_exist<int> ("label");
  }

  // Describes how to process one line (= one point object)
  void process_line (CGAL::Ply_reader& reader)
  {
    double x = 0., y = 0., z = 0.;
    int l = 0;

    reader.assign (x, "x");
    reader.assign (y, "y");
    reader.assign (z, "z");
    reader.assign (l, "label");

    points.push_back (Point (x, y, z));
    labels.push_back(std::size_t(l));
  }

};


int main (int argc, char** argv)
{
  std::string filename (argc > 1 ? argv[1] : "data/b9_training.ply");
  std::ifstream in (filename.c_str());
  std::vector<Point> pts;
  std::vector<std::size_t> ground_truth;

  std::cerr << "Reading input" << std::endl;
  My_ply_interpreter interpreter (pts, ground_truth);
  if (!in
      || !(CGAL::read_ply_custom_points (in, interpreter, Kernel())))
  {
    std::cerr << "Error: cannot read " << filename << std::endl;
    return EXIT_FAILURE;
  }

  Feature_set features;
  
  std::cerr << "Generating features" << std::endl;
  CGAL::Real_timer t;
  t.start();
  Feature_generator generator (features, pts, Pmap(),
                               5);  // using 5 scales
  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;

  // Add types
  Label_set labels;
  Label_handle ground = labels.add ("ground");
  Label_handle vege = labels.add ("vegetation");
  Label_handle roof = labels.add ("roof");
  Label_handle facade = labels.add ("facade");

  Classification_predicate predicate (labels, features);
  
  std::cerr << "Training" << std::endl;
  t.reset();
  t.start();
  predicate.train (ground_truth);
  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;

  t.reset();
  t.start();
  std::vector<std::size_t> label_indices;
  Classif::classify_with_graphcut<CGAL::Sequential_tag>
    (pts, Pmap(), labels, predicate,
     generator.neighborhood().k_neighbor_query(12),
     0.2, 10, label_indices);
  t.stop();
  std::cerr << "Classification with graphcut done in " << t.time() << " second(s)" << std::endl;

  std::cerr << "Precision, recall, F1 scores and IoU:" << std::endl;
  Classif::Evaluation eval (labels, ground_truth, label_indices);
  
  for (std::size_t i = 0; i < labels.size(); ++ i)
  {
    std::cerr << " * " << labels[i]->name() << ": "
              << eval.precision(labels[i]) << " ; "
              << eval.recall(labels[i]) << " ; "
              << eval.f1_score(labels[i]) << " ; "
              << eval.intersection_over_union(labels[i]) << std::endl;
  }

  std::cerr << "Accuracy = " << eval.accuracy() << std::endl
            << "Mean F1 score = " << eval.mean_f1_score() << std::endl
            << "Mean IoU = " << eval.mean_intersection_over_union() << std::endl;
  
  std::ofstream f ("classification.ply");
  f.precision(18);
  f << "ply" << std::endl
    << "format ascii 1.0" << std::endl
    << "element vertex " << pts.size() << std::endl
    << "property float x" << std::endl
    << "property float y" << std::endl
    << "property float z" << std::endl
    << "property uchar red" << std::endl
    << "property uchar green" << std::endl
    << "property uchar blue" << std::endl
    << "end_header" << std::endl;
  
  for (std::size_t i = 0; i < pts.size(); ++ i)
    {
      f << pts[i] << " ";
      
      Label_handle label = labels[label_indices[i]];
      if (label == ground)
        f << "245 180 0" << std::endl;
      else if (label == vege)
        f << "0 255 27" << std::endl;
      else if (label == roof)
        f << "255 0 170" << std::endl;
      else if (label == facade)
        f << "128 128 128" << std::endl;
      else
        {
          f << "0 0 0" << std::endl;
          std::cerr << "Error: unknown classification label" << std::endl;
        }
    }

  std::cerr << "All done" << std::endl;
  
  return EXIT_SUCCESS;
}
