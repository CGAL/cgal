#define CGAL_DELAUNAY_3_USE_SUBDETERMINANTS 1
// #define CGAL_PROFILE
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Timer.h>
#include <vector>
#include <array>
#include <iostream>
#include <string>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_3<K>                        Vb;
#ifdef CGAL_DELAUNAY_3_USE_SUBDETERMINANTS
typedef CGAL::Triangulation_cell_base_with_info_3<std::array<double,4>, K>    Cbb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K,Cbb>             Cb;
#else
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
#endif
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                      Delaunay;
typedef Delaunay::Point                                             Point;
typedef CGAL::Timer                                                 Timer;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/ocean_r.xyz");
  std::ifstream in(filename.c_str());
  std::vector<Point> points;
  Point p;

  while(in >> p){
    points.push_back(p);
  }

  std::cout << points.size() << " points read\n";

  Timer timer;
  timer.start();
  {
    Delaunay dt;
    dt.insert(points.begin(), points.end());
  }
  std::cerr << timer.time() << " sec" << std::endl;
  std::cout << "done" << std::endl;
  return 0;
}

