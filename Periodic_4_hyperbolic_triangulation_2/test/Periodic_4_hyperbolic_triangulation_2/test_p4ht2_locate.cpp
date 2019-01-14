#include <CGAL/Cartesian.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/internal/Periodic_4_hyperbolic_triangulation_dummy_14.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <iostream>

typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<>               Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>                Triangulation;
typedef Triangulation::Face_handle                                                  Face_handle;
typedef Triangulation::Vertex_handle                                                Vertex_handle;
typedef Triangulation::Locate_type                                                  Locate_type;
typedef Triangulation::Hyperbolic_translation                                       Hyperbolic_translation;
typedef Triangulation::Point                                                        Point;

std::ostream& operator<<(std::ostream& s, const Locate_type& lt)
{
  switch(lt)
  {
    case Triangulation::VERTEX:
      s << "VERTEX";
      break;
    case Triangulation::FACE:
      s << "FACE";
      break;
    case Triangulation::EDGE:
      s << "EDGE";
      break;
  }

  return s;
}

int main(int, char**)
{
  Triangulation tr;

  assert(tr.is_valid());

  Locate_type lt;
  int li;
  Face_handle fh;

  std::cout << "---- locating dummy points (all should be vertices) ----" << std::endl;
  for(int j=0; j<14; ++j) {
    Point query = tr.get_dummy_point(j);
    fh = tr.hyperbolic_locate(query, lt, li);
    assert(lt == Triangulation::VERTEX);
    std::cout << "   dummy point " << j << ": OK " << std::endl;
  }

  std::cout << "---- locating the midpoint of a Euclidean segment ----" << std::endl;
  Point p1 = tr.get_dummy_point(0), p2 = tr.get_dummy_point(1);
  Point query = midpoint(p1, p2);
  fh = tr.hyperbolic_locate(query, lt, li);
  assert(lt == Triangulation::EDGE);
  std::cout << "   located as edge OK" << std::endl;

  std::cout << "---- inserting a single point and locating it ----" << std::endl;
  Vertex_handle v = tr.insert(Point(-0.4, -0.1));
  fh = tr.hyperbolic_locate(v->point(), lt, li);
  assert(lt == Triangulation::VERTEX);
  std::cout << "   located as vertex OK" << std::endl;

  // TODO: add a test case for a circular edge!

  return EXIT_SUCCESS;
}
