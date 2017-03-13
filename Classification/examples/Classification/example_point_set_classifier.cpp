#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//#define CGAL_CLASSIFICATION_VERBOSE

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_classifier.h>
#include <CGAL/Classification/Trainer.h>
#include <CGAL/IO/read_ply_points.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point> Point_range;
typedef CGAL::Identity_property_map<Point> Pmap;

typedef CGAL::Point_set_classifier<Kernel, Point_range, Pmap> Point_set_classifier;
typedef CGAL::Classification::Trainer<Point_range, Pmap> Trainer;

/*
 This interpreter is used to read a PLY input that contains training
 attributes (with the PLY "label" property).
*/
class My_ply_interpreter
{
  std::vector<Point>& points;
  std::vector<int>& labels;
    
public:
  My_ply_interpreter (std::vector<Point>& points,
                      std::vector<int>& labels)
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
    labels.push_back(l);
  }

};


int main (int argc, char** argv)
{
  std::string filename (argc > 1 ? argv[1] : "data/b9_training.ply");
  std::ifstream in (filename.c_str());
  std::vector<Point> pts;
  std::vector<int> labels;

  std::cerr << "Reading input" << std::endl;
  My_ply_interpreter interpreter (pts, labels);
  if (!in
      || !(CGAL::read_ply_custom_points (in, interpreter, Kernel())))
    {
      std::cerr << "Error: cannot read " << filename << std::endl;
      return EXIT_FAILURE;
    }

  Point_set_classifier psc (pts, Pmap());
  
  std::cerr << "Generating features" << std::endl;
  CGAL::Real_timer t;
  t.start();
  psc.generate_features (5); // Using 5 scales
  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;
  
  // Add types to PSC
  CGAL::Classification::Label_handle ground
    = psc.add_label ("ground");
  CGAL::Classification::Label_handle vege
    = psc.add_label ("vegetation");
  CGAL::Classification::Label_handle roof
    = psc.add_label ("roof");
  CGAL::Classification::Label_handle facade
    = psc.add_label ("facade");
  
  Trainer trainer (psc);

  // Set training sets
  std::size_t nb_inliers = 0;
  for (std::size_t i = 0; i < labels.size(); ++ i)
    {
      switch (labels[i])
        {
        case 0:
          trainer.set_inlier(vege, i);
          ++ nb_inliers;
          break;
        case 1:
          trainer.set_inlier(ground, i);
          ++ nb_inliers;
          break;
        case 2:
          trainer.set_inlier(roof, i);
          ++ nb_inliers;
          break;
        case 3:
          trainer.set_inlier(facade, i);
          ++ nb_inliers;
          break;
        default:
          break;
        }
    }

  std::cerr << "Training using " << nb_inliers << " inliers" << std::endl;
  t.reset();
  t.start();
  trainer.train (400); // 800 trials
  t.stop();
  std::cerr << "Done in " << t.time() << " second(s)" << std::endl;

  std::cerr << "Precision, recall, F1 scores and IoU:" << std::endl;
  for (std::size_t i = 0; i < psc.number_of_labels(); ++ i)
    {
      std::cerr << " * " << psc.label(i)->name() << ": "
                << trainer.precision(psc.label(i)) << " ; "
                << trainer.recall(psc.label(i)) << " ; "
                << trainer.f1_score(psc.label(i)) << " ; "
                << trainer.intersection_over_union(psc.label(i)) << std::endl;
    }

  std::cerr << "Accuracy = " << trainer.accuracy() << std::endl
            << "Mean F1 score = " << trainer.mean_f1_score() << std::endl
            << "Mean IoU = " << trainer.mean_intersection_over_union() << std::endl;
  
  t.reset();
  t.start();
  psc.run_with_graphcut (psc.neighborhood().k_neighbor_query(12), 0.5);
  t.stop();
  std::cerr << "One graphcut done in " << t.time() << " second(s)" << std::endl;
  
  // Save the output in a colored PLY format
  {
    std::ofstream f ("classification_one.ply");
    f.precision(18);
    psc.write_classification_to_ply (f);
  }

  t.reset();
  t.start();
  psc.run_with_graphcut (psc.neighborhood().k_neighbor_query(12), 0.5, 30);
  t.stop();
  std::cerr << std::size_t(pts.size() / 25000) << " graphcuts done in " << t.time() << " second(s)" << std::endl;
  
  // Save the output in a colored PLY format
  {
    std::ofstream f ("classification_several.ply");
    f.precision(18);
    psc.write_classification_to_ply (f);
  }
  
  /// Save the configuration to be able to reload it later
  std::ofstream fconfig ("config.xml");
  psc.save_configuration (fconfig);
  
  std::cerr << "All done" << std::endl;
  
  return EXIT_SUCCESS;
}
