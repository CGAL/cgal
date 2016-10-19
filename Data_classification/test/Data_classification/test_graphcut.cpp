#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/Data_classification/Helper.h>
#include <CGAL/IO/read_ply_points.h>

#define USE_POINT_SET_3

#ifdef USE_POINT_SET_3
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#endif

//typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

#ifdef USE_POINT_SET_3
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::iterator Iterator;
typedef Point_set::Point_map Pmap;
#else
typedef std::vector<Point>::iterator Iterator;
typedef CGAL::Identity_property_map<Point> Pmap;
#endif

typedef CGAL::Point_set_classification<Kernel, Iterator, Pmap> Classification;
typedef CGAL::Data_classification::Helper<Kernel, Iterator, Pmap> Helper;

int main (int argc, char** argv)
{

  std::string filename (argc > 1 ? argv[1] : "data/portland_mini.ply");
  std::ifstream in (filename.c_str());

#ifdef USE_POINT_SET_3
  Point_set pts;
#else
  std::vector<Point> pts;
#endif

  std::cerr << "[READING INPUT]" << std::endl;
  if (!in
#ifdef USE_POINT_SET_3
      || !(CGAL::read_ply_point_set (in, pts)))
#else
      || !(CGAL::read_ply_points (in, std::back_inserter (pts))))
#endif
    {
      std::cerr << "Error: cannot read " << filename << std::endl;
      return EXIT_FAILURE;
    }

  std::string configname (argc > 2 ? argv[2] : "data/portland_mini.xml");

  std::cerr << "[LOADING CONFIG]" << std::endl;
#ifdef USE_POINT_SET_3

  Classification psc (pts.begin (), pts.end(), pts.point_map());

  Helper helper (configname.c_str(), // load XML config
                 psc,
                 pts.begin(), pts.end(), pts.point_map());
#else
  Classification psc (pts.begin (), pts.end(), Pmap());

  Helper helper (configname.c_str(), // load XML config
                 psc,
                 pts.begin(), pts.end(), Pmap());
#endif
  std::cerr << "[RUNNING WITH GRAPHCUT]" << std::endl;

  psc.run_with_graphcut (helper.neighborhood(), 0.5);

  std::cerr << std::endl << "[WRITING OUTPUT]" << std::endl;

  std::ofstream f ("output.ply");
  f.precision(18);
  srand (time (NULL)); // to get different colors everytime
#ifdef USE_POINT_SET_3
  helper.write_ply (f, pts.begin(), pts.end(), pts.point_map(), psc);
#else
  helper.write_ply (f, pts.begin(), pts.end(), Pmap(), psc);
#endif

  std::cerr << "[ALL DONE]" << std::endl;

  return EXIT_SUCCESS;
}
