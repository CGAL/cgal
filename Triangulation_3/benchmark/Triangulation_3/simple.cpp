//#define CGAL_PROFILE
//#define CGAL_USE_SSE2_FABS
//#define CGAL_USE_SSE2_MAX
//#define CGAL_MSVC_USE_STD_FABS  // use this one with precise

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include <CGAL/Delaunay_triangulation_3.h>


#include <CGAL/Timer.h>
#include <iostream>
#include <string>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Delaunay_triangulation_3<K>                    DT;
typedef DT::Point                                            Point_3;
typedef CGAL::Timer                                          Timer;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/ocean_r.xyz");
  std::ifstream in(filename.c_str());
  std::vector<Point_3> points;
  Point_3 p, q;

  while(in >> p ){
    points.push_back(p);
  }

  std::cout << points.size() << " points read\n";

  Timer timer;
  timer.start();
  size_t N = 0;
  for(int i = 0; i < 1; i++){
    DT dt;
    dt.insert(points.begin(), points.end());
    N += dt.number_of_cells();
  }
  timer.stop();

  std::cerr << N << std::endl << timer.time() << " sec" << std::endl;
  return 0;
}
