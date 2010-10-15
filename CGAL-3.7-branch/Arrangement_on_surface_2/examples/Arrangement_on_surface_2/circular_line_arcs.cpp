//! \file examples/Arrangement_on_surface_2/circular_line_arc.cpp
// Using the circular line arc traits.

#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Arr_circular_line_arc_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <boost/variant.hpp>

#include <CGAL/Random.h>


typedef CGAL::Quotient<CGAL::MP_Float>                      NT;
typedef CGAL::Cartesian<NT>                                 Linear_k;
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT>          Algebraic_k;
typedef CGAL::Circular_kernel_2<Linear_k,Algebraic_k>       Circular_k;

typedef Circular_k::Point_2                                 Point_2;
typedef Circular_k::Circle_2                                Circle_2;
typedef Circular_k::Circular_arc_2                          Circular_arc_2;
typedef Circular_k::Line_arc_2                              Line_arc_2;

typedef boost::variant< Circular_arc_2, Line_arc_2>         Arc_2;
typedef std::vector< Arc_2>                                 ArcContainer;

typedef CGAL::Arr_circular_line_arc_traits_2<Circular_k>    Traits;

typedef CGAL::Arrangement_2<Traits>                         Arrangement;
typedef CGAL::Arr_naive_point_location<Arrangement>         Point_location;

int main()
{
  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  std::cout << "random_seed = " << random_seed << std::endl;
  CGAL::Random theRandom(random_seed);
  int random_max = 128;
  int random_min = -128;
  ArcContainer ac;
  int x1, y1, x2, y2;

  for (int i = 0; i < 10; i++) {
    x1 = theRandom.get_int(random_min,random_max);
    y1 = theRandom.get_int(random_min,random_max);
    do{
      x2 = theRandom.get_int(random_min,random_max);
      y2 = theRandom.get_int(random_min,random_max);
    } while((x1 == x2) && (y1 == y2));

    std::cout << x1 << " " << y1 << " " << x2 << " " << y2 << std::endl;
    boost::variant< Circular_arc_2, Line_arc_2 > v =
      Line_arc_2(Point_2(x1,y1), Point_2(x2,y2));
    ac.push_back( v);
  }

  for (int i = 0; i < 10; i++) {
    x1 = theRandom.get_int(random_min,random_max);
    y1 = theRandom.get_int(random_min,random_max);
    boost::variant< Circular_arc_2, Line_arc_2 > v =
      Circle_2( Point_2(x1,y1), x1*x1 + y1*y1);
    ac.push_back(v);
  }

  Arrangement arr;
  Point_location _pl(arr);
  for (ArcContainer::const_iterator it = ac.begin(); it != ac.end(); ++it) {
    //insert(arr,_pl,*it);
    insert(arr, *it, _pl);
  };

  return 0;
}
