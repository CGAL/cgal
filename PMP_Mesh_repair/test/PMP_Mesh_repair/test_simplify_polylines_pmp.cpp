#include <CGAL/Polygon_mesh_processing/internal/simplify_polyline.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <vector>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;

namespace PMP_exp=CGAL::Polygon_mesh_processing::experimental;

int main(int argc, char **argv)
{
  for (int i=1; i<argc; i+=2)
  {
    if (i+1>=argc)
    {
      std::cerr << "Usage: " << argv[0] << " in.polylines_0.txt max_squared_frechet_distance_0 ... in.polylines_n.txt max_squared_frechet_distance_n \n";
      return 1;
    }

    std::ifstream in(argv[i]);
    const double max_squared_frechet_distance = atof(argv[i+1]);

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

    std::cout << "Polyline of " << n << " points simplified (ITERATIVE) into a polyline of " << simplified_polyline.size() << " points\n";

    std::ofstream out("out_iterative_"+std::to_string((i-1)/2)+".polylines.txt");
    out << simplified_polyline.size();
    for (auto p : simplified_polyline)
      out << " " << p;
    out << "\n";
    out.close();

    // second run
    simplified_polyline.clear();
    PMP_exp::simplify_polyline(polyline, simplified_polyline, max_squared_frechet_distance,
                               CGAL::parameters::algorithm(PMP_exp::DOUGLAS_PEUCKER));

    std::cout << "Polyline of " << n << " points simplified (DOUGLAS_PEUCKER) into a polyline of " << simplified_polyline.size() << " points\n";

    out.open("out_dp_"+std::to_string((i-1)/2)+".polylines.txt");
    out << simplified_polyline.size();
    for (auto p : simplified_polyline)
      out << " " << p;
    out << "\n";
  }

  return 0;
}
