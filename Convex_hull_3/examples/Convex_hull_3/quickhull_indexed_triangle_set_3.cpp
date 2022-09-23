#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>

#include <vector>
#include <array>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::Point_3                                Point_3;


int main(int argc, char* argv[])
{
  std::ifstream in( (argc>1)? argv[1] : "data/cube.xyz");
  std::vector<Point_3> points;

  std::vector<Point_3> vertices;
  std::vector<std::array<int,3> > faces;

  Point_3 p;
  while(in >> p){
    points.push_back(p);
  }


  CGAL::convex_hull_3(points.begin(), points.end(), vertices, faces);

  std::cout << vertices.size() << "  " << faces.size() << std::endl;

  return 0;
}
