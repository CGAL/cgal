#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_classifier.h>
#include <CGAL/Classification/Trainer.h>
#include <CGAL/IO/read_ply_points.h>

#include <CGAL/Timer.h>

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
  
  std::cerr << "Generating attributes" << std::endl;
  psc.generate_attributes (5); // Using 5 scales
  
  // Add types to PSC
  CGAL::Classification::Type_handle ground
    = psc.add_classification_type ("ground");
  CGAL::Classification::Type_handle vege
    = psc.add_classification_type ("vegetation");
  CGAL::Classification::Type_handle roof
    = psc.add_classification_type ("roof");
  
  Trainer trainer (psc);

  // Set training sets
  for (std::size_t i = 0; i < labels.size(); ++ i)
    {
      switch (labels[i])
        {
        case 0:
          trainer.set_inlier(ground, i);
          break;
        case 1:
          trainer.set_inlier(vege, i);
          break;
        case 2:
          trainer.set_inlier(roof, i);
          break;
        default:
          break;
        }
    }

  std::cerr << "Training" << std::endl;
  trainer.train (800); // 800 trials
  
  psc.run_with_graphcut (psc.neighborhood().k_neighbor_query(12), 0.5);
  
  // Save the output in a colored PLY format
  std::ofstream f ("classification.ply");
  f.precision(18);

  psc.write_classification_to_ply (f);

  /// Save the configuration to be able to reload it later
  std::ofstream fconfig ("config.xml");
  psc.save_configuration (fconfig);
  
  std::cerr << "All done" << std::endl;
  
  return EXIT_SUCCESS;
}
