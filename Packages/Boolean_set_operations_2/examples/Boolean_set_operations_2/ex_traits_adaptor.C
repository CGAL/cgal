//! \file examples/Boolean_set_operations_2/ex_traits_adaptor.C
// Using the traits adaptor to generate a traits of polylines.

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Gps_traits_adaptor_2.h>

#include <list>

typedef CGAL::Gmpq                                      NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Arr_segment_traits;
typedef CGAL::Arr_polyline_traits_2<Arr_segment_traits> Arr_traits;
typedef CGAL::General_polygon_2<Arr_traits>             Polygon;
typedef CGAL::Gps_traits_adaptor_2<Arr_traits,Polygon>  Traits;
typedef CGAL::General_polygon_set_2<Traits>             Polygon_set;
typedef Traits::Polygon_with_holes_2                    Polygon_with_holes;
typedef Traits::Curve_2                                 Curve;
typedef Traits::X_monotone_curve_2                      X_monotone_curve;
typedef Traits::Point_2                                 Point;

void polyline_2_polygon(Curve polyline, Polygon & polygon)
{
  Traits traits;
  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object()(polyline, std::back_inserter(objects));
  std::list<CGAL::Object>::iterator i;
  for (i = objects.begin(); i != objects.begin(); ++i) {
    X_monotone_curve xcurve;
    CGAL::assign(xcurve, *i);
    polygon.push_back(xcurve);
  }
}

int main(int argc, char * argv[])
{
  Polygon pn1, pn2;
  Curve pl1, pl2;

  pl1.push_back(Point(0,0));     pl1.push_back(Point(1.5,1.5));
  pl1.push_back(Point(2.5,0.5)); pl1.push_back(Point(3.5,1.5));
  pl1.push_back(Point(5,0));
  pl2.push_back(Point(0,2));     pl2.push_back(Point(5,2));
  pl2.push_back(Point(3.5,0.5)); pl2.push_back(Point(2.5,1.5));
  pl2.push_back(Point(1.5,0.5));
  
  polyline_2_polygon(pl1, pn1);
  polyline_2_polygon(pl2, pn2);
  std::list<Polygon_with_holes> result;
  CGAL::intersection(pn1, pn2, std::back_inserter(result));

  std::copy(result.begin(), result.end(),       // export to standard output
            std::ostream_iterator<Polygon_with_holes>(std::cout, "\n"));
  std::cout << std::endl;
  
  return 0;
}
