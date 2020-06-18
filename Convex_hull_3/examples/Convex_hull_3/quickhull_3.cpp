#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Point_3                                Point_3;
typedef CGAL::Surface_mesh<Point_3>               Surface_mesh;


int main(int argc, char* argv[])
{
  std::ifstream in( (argc>1)? argv[1] : "data/cube.xyz");
  std::vector<Point_3> points;
  Point_3 p;
  while(in >> p){
    points.push_back(p);
  }

  // define polyhedron to hold convex hull
  Polyhedron_3 poly;

  // compute convex hull of non-collinear points
  CGAL::convex_hull_3(points.begin(), points.end(), poly);

  std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices" << std::endl;

  Surface_mesh sm;
  CGAL::convex_hull_3(points.begin(), points.end(), sm);

  std::cout << "The convex hull contains " << num_vertices(sm) << " vertices" << std::endl;

  return 0;
}
