#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Indexed_triangle_set.h>
#include <CGAL/convex_hull_3.h>

#include <vector>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::Point_3                                Point_3;
typedef CGAL::Indexed_triangle_set<Point_3>               Surface_mesh;


int main(int argc, char* argv[])
{
  std::ifstream in( (argc>1)? argv[1] : "data/cube.xyz");
  std::vector<Point_3> points;
  Point_3 p;
  while(in >> p){
    points.push_back(p);
  }

  Surface_mesh sm;
  CGAL::convex_hull_3(points.begin(), points.end(), sm);

  std::cout << sm << std::endl;
  //std::cout << "The convex hull contains " << num_vertices(sm) << " vertices" << std::endl;

  return 0;
}
