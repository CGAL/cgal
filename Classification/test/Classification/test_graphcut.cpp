#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Classifier.h>
#include <CGAL/Classification/Helper.h>
#include <CGAL/IO/read_ply_points.h>

#include <CGAL/Timer.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

//typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::iterator Iterator;
typedef Point_set::Point_map Pmap;

typedef CGAL::Classifier<Iterator, Pmap> Classification;
typedef CGAL::Classification::Helper<Kernel, Iterator, Pmap> Helper;

int main (int argc, char** argv)
{

  std::string filename (argc > 1 ? argv[1] : "data/portland_mini.ply");
  std::ifstream in (filename.c_str());

  Point_set pts;
  std::cerr << "[READING INPUT]" << std::endl;
  if (!in
      || !(CGAL::read_ply_point_set (in, pts)))
    {
      std::cerr << "Error: cannot read " << filename << std::endl;
      return EXIT_FAILURE;
    }

  std::string configname (argc > 2 ? argv[2] : "data/portland_mini.xml");

  std::cerr << "[LOADING CONFIG]" << std::endl;

  Classification psc (pts.begin (), pts.end(), pts.point_map());

  Helper helper (configname.c_str(), // load XML config
                 psc,
                 pts.begin(), pts.end(), pts.point_map());

  std::cerr << "[RUNNING WITH GRAPHCUT]" << std::endl;

  CGAL::Timer t;
  t.start();
  psc.run_with_graphcut (helper.neighborhood(), 0.5);
  t.stop();
  std::cerr << "Grapcut computed in " << t.time() << " second(s)" << std::endl;

  std::cerr << std::endl << "[WRITING OUTPUT]" << std::endl;

  std::ofstream f ("output.ply");
  f.precision(18);
  srand (time (NULL)); // to get different colors everytime

  helper.write_ply (f, pts.begin(), pts.end(), pts.point_map(), psc);

  std::cerr << "[ALL DONE]" << std::endl;

  return EXIT_SUCCESS;
}
