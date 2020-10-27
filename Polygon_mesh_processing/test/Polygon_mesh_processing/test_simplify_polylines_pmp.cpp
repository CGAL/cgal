#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/internal/simplify_polyline.h>

#include <vector>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;

namespace PMP_exp=CGAL::Polygon_mesh_processing::experimental;

int main(int argc, char **argv)
{
  if (argc!=3)
  {
    std::cerr << "Usage: " << argv[0] << " in.polylines.txt max_squared_frechet_distance\n";
    return 1;
  }

  std::ifstream in(argv[1]);
  const double max_squared_frechet_distance = atof(argv[2]);

  //TODO read more than 1 poly or use CGAL I/O function
  int n;
  in >> n;

  std::vector<Point_3> polyline(n);
  for(int i=0; i<n; ++i)
    in >> polyline[i];

  // first run
  std::vector<Point_3> simplified_polyline;
  PMP_exp::simplify_polyline(polyline, simplified_polyline, max_squared_frechet_distance, 
                             CGAL::parameters::algorithm(PMP_exp::ITERATIVE));

  std::cout << "Polyline of " << n << " points simplified into a polyline of " << simplified_polyline.size() << " points\n";

  std::ofstream out("out_iterative.polylines.txt");
  out << simplified_polyline.size();
  for (auto p : simplified_polyline)
    out << " " << p;
  out << "\n";
  out.close();

  // second run
  simplified_polyline.clear();
  PMP_exp::simplify_polyline(polyline, simplified_polyline, max_squared_frechet_distance, 
                             CGAL::parameters::algorithm(PMP_exp::DOUGLAS_PEUCKER));

  std::cout << "Polyline of " << n << " points simplified into a polyline of " << simplified_polyline.size() << " points\n";

  out.open("out_dp.polylines.txt");
  out << simplified_polyline.size();
  for (auto p : simplified_polyline)
    out << " " << p;
  out << "\n";

  return 0;
}
