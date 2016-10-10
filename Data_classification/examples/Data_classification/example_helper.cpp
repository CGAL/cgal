#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/Data_classification/Helper.h>
#include <CGAL/IO/read_ply_points.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point>::iterator Iterator;
typedef CGAL::Identity_property_map<Point> Pmap;

typedef CGAL::Point_set_classification<Kernel, Iterator, Pmap> Classification;
typedef CGAL::Data_classification::Helper<Kernel, Iterator, Pmap> Helper;

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
  std::string filename (argc > 1 ? argv[1] : "data/b9.ply");
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

  double grid_resolution = 1.;
  double radius_neighbors = 5.;
  double radius_dtm = 30;

  Classification psc (pts.begin (), pts.end(), Pmap());
  
  std::cerr << "Using helper " << std::endl;

  // for (std::size_t i = 0; i < 3; ++ i)
  //   {
      Helper helper (pts.begin(), pts.end(), Pmap(),
                     grid_resolution, radius_neighbors, radius_dtm);

      helper.generate_attributes (psc, pts.begin(), pts.end(), Pmap());
      //    }
  
  // Add types to PSC
  CGAL::Data_classification::Type_handle ground
    = psc.add_classification_type ("ground");
  CGAL::Data_classification::Type_handle vege
    = psc.add_classification_type ("vegetation");
  CGAL::Data_classification::Type_handle roof
    = psc.add_classification_type ("roof");
  CGAL::Data_classification::Type_handle facade
    = psc.add_classification_type ("facade");

  for (std::size_t i = 0; i < labels.size(); ++ i)
    {
      switch (labels[i])
        {
        case 0:
          psc.add_training_index(vege, i);
          break;
        case 1:
          psc.add_training_index(ground, i);
          break;
        case 2:
          psc.add_training_index(roof, i);
          break;
        case 3:
          psc.add_training_index(facade, i);
          break;
        default:
          break;
        }
    }


  //  psc.set_multiplicative(true);
  std::cerr << "Training" << std::endl;
  psc.training(1000);

  //  psc.run_with_graphcut (helper.neighborhood(), 0.5);
  psc.run();
  
  // Save the output in a colored PLY format

  std::ofstream f ("classification.ply");
  f.precision(18);
  f << "ply" << std::endl
    << "format ascii 1.0" << std::endl
    << "element vertex " << pts.size() << std::endl
    << "property double x" << std::endl
    << "property double y" << std::endl
    << "property double z" << std::endl
    << "property uchar red" << std::endl
    << "property uchar green" << std::endl
    << "property uchar blue" << std::endl
    << "property int label" << std::endl
    << "end_header" << std::endl;
  
  for (std::size_t i = 0; i < pts.size(); ++ i)
    {
      f << pts[i] << " ";
      
      CGAL::Data_classification::Type_handle type = psc.classification_type_of (i);
      if (type == vege)
        f << "0 255 27 0" << std::endl;
      else if (type == ground)
        f << "245 180 0 1" << std::endl;
      else if (type == roof)
        f << "255 0 170 2" << std::endl;
      else if (type == facade)
        f << "100 0 255 3" << std::endl;
      else
        {
          f << "0 0 0 -1" << std::endl;
          //          std::cerr << "Error: unknown classification type" << std::endl;
        }
    }
  helper.save ("config.xml", psc);
  
  std::cerr << "All done" << std::endl;
  
  return EXIT_SUCCESS;
}
