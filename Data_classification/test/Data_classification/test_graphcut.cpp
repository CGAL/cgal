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
typedef std::vector<Point>::iterator Iterator;
typedef CGAL::Identity_property_map<Point> Pmap;

typedef CGAL::Point_set_classification<Kernel, Iterator, Pmap> Classification;
typedef CGAL::Data_classification::Helper<Kernel, Iterator, Pmap> Helper;

int main (int argc, char** argv)
{

  std::string filename (argc > 1 ? argv[1] : "data/portland_mini.ply");
  std::ifstream in (filename.c_str());
  std::vector<Point> pts;

  std::cerr << "[READING INPUT]" << std::endl;
  if (!in
      || !(CGAL::read_ply_points (in, std::back_inserter (pts))))
    {
      std::cerr << "Error: cannot read " << filename << std::endl;
      return EXIT_FAILURE;
    }

  std::string configname (argc > 2 ? argv[2] : "data/portland_mini.xml");

  std::cerr << "[LOADING CONFIG]" << std::endl;

  Classification psc (pts.begin (), pts.end(), Pmap());

  Helper helper (configname.c_str(), // load XML config
                 psc,
                 pts.begin(), pts.end(), Pmap());

  std::cerr << "[RUNNING WITH GRAPHCUT]" << std::endl;

  psc.run_with_graphcut (helper.neighborhood(), 0.5);

  std::cerr << std::endl << "[WRITING OUTPUT]" << std::endl;

  std::ofstream f ("output.ply");
  f.precision(18);
  srand (time (NULL)); // to get different colors everytime
  helper.write_ply (f, pts.begin(), pts.end(), Pmap(), psc);

  std::cerr << "[ALL DONE]" << std::endl;

  return EXIT_SUCCESS;
}
