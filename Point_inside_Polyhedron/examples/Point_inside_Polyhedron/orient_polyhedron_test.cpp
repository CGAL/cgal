#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Orient_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Timer.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polygon_soup_orienter<Polyhedron> Polygon_soup_orienter;

int main() {
  std::ifstream input("elephant_soup.off");
  Polygon_soup_orienter soup;
  if ( !input || !(input >> soup)){
    std::cerr << "Error: can not read file.";
    return 1;
  }

  Polyhedron poly;
  bool oriented = soup.orient(poly);

  std::cerr << (oriented ? "Oriented." : "Not oriented") << std::endl;
  if(oriented) {
    std::ofstream out("elephant_oriented.off");
    out << poly;
    out.close();
  }
}