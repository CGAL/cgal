#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Point_3                                Point_3;
typedef K::Segment_3                              Segment_3;
typedef K::Triangle_3                             Triangle_3;


int main(int argc, char* argv[])
{
  std::ifstream in( (argc>1)? argv[1] : "data/cube.xyz");
  std::vector<Point_3> points;
  Point_3 p;

  while(in >> p){
    points.push_back(p);
  }

  CGAL::Object obj;

  // compute convex hull of non-collinear points
  CGAL::convex_hull_3(points.begin(), points.end(), obj);

  if(const Point_3* p = CGAL::object_cast<Point_3>(&obj)){
    std::cout << "Point " << *p << std::endl;
  }
  else if(const Segment_3* s = CGAL::object_cast<Segment_3>(&obj)){
    std::cout << "Segment " << *s << std::endl;
  }
  else if(const Triangle_3* t = CGAL::object_cast<Triangle_3>(&obj)){
    std::cout << "Triangle " << *t << std::endl;
  }
  else  if(const Polyhedron_3* poly = CGAL::object_cast<Polyhedron_3>(&obj)){
    std::cout << "Polyhedron\n " << *poly << std::endl;
  std::cout << "The convex hull contains " << poly->size_of_vertices() << " vertices" << std::endl;

  }
  else {
    std::cout << "something else"<< std::endl;
  }
  return 0;
}
