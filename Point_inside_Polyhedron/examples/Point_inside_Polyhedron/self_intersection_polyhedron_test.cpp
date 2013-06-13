#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Self_intersection_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Timer.h>
#include <boost/lexical_cast.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef CGAL::Polyhedron_3<K> Polyhedron;

int main() {
  std::ifstream input("camel.off");
  Polyhedron poly;
  if ( !input || !(input >> poly) || poly.empty() ){
    std::cerr << "Error: can not read file.";
    return 1;
  }
  
  CGAL::Timer timer;
  timer.start();

  std::vector<Triangle> intersected_tris;
  typedef std::back_insert_iterator<std::vector<Triangle> > OutputIterator;
  CGAL::self_intersect<Polyhedron, K, OutputIterator>(poly, OutputIterator(intersected_tris));

  std::cerr << "Self-intersection test took " << timer.time() << " sec." << std::endl;
  std::cerr << intersected_tris.size() << " triangles are intersecting." << std::endl;

  timer.reset();
  bool intersecting = CGAL::self_intersect<Polyhedron, K>(poly);

  std::cerr << "Is self-intersection test took " << timer.time() << " sec." << std::endl;
  std::cerr << (intersecting ? "There is an self-intersection." : "There is no self-intersection.") << std::endl;
}