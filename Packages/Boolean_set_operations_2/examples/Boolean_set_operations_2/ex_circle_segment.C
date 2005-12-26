//! \file examples/Boolean_set_operations_2/ex_circle_segment.C
// Handling circles and linear segments concurrently.

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/General_polygon_with_holes_2.h>

#include <list>

typedef CGAL::Gmpq                                      NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Gps_circle_segment_traits_2<Kernel>       Traits;
typedef CGAL::General_polygon_set_2<Traits>             Polygon_set_2;
typedef Traits::Polygon_2                               Polygon_2;
typedef Traits::Polygon_with_holes_2                    Polygon_with_holes_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef Kernel::Point_2                                 Point_2;
typedef Kernel::Circle_2                                Circle_2;

void circle_2_polygon(Circle_2 circle, Polygon_2 & polygon)
{
  Traits traits;
  Curve_2 curve(circle);
  std::list<CGAL::Object> objects;
  traits.make_x_monotone_2_object()(curve, std::back_inserter(objects));
  std::list<CGAL::Object>::iterator i = objects.begin();
  X_monotone_curve_2 xcurve;
  CGAL::assign(xcurve, *i++);   polygon.push_back(xcurve);
  CGAL::assign(xcurve, *i);     polygon.push_back(xcurve);
}

void rect_2_polygon(Point_2 p1, Point_2 p2, Point_2 p3, Point_2 p4,
                    Polygon_2 & polygon)
{
  X_monotone_curve_2 s1(p1, p2);  polygon.push_back(s1);
  X_monotone_curve_2 s2(p2, p3);  polygon.push_back(s2);
  X_monotone_curve_2 s3(p3, p4);  polygon.push_back(s3);
  X_monotone_curve_2 s4(p4, p1);  polygon.push_back(s4);
}

int main(int argc, char * argv[])
{
  Polygon_2 pc1, pc2, pc3, pc4, pr1, pr2, pr3, pr4;
  Polygon_set_2 ps;

  // Input non-intersecting circles:
  circle_2_polygon(Circle_2(Point_2(1, 1), 1), pc1);  ps.insert(pc1);
  circle_2_polygon(Circle_2(Point_2(5, 1), 1), pc2);  ps.insert(pc2);
  circle_2_polygon(Circle_2(Point_2(5, 5), 1), pc3);  ps.insert(pc3);
  circle_2_polygon(Circle_2(Point_2(1, 5), 1), pc4);  ps.insert(pc4);
  
  // Perform union with rectangles incrementally:
  rect_2_polygon(Point_2(1,0), Point_2(5,0), Point_2(5,2), Point_2(1,2), pr1);
  ps.join(pr1);

  rect_2_polygon(Point_2(1,4), Point_2(5,4), Point_2(5,6), Point_2(1,6), pr2);
  ps.join(pr2);

  rect_2_polygon(Point_2(0,1), Point_2(0,5), Point_2(2,5), Point_2(2,1), pr3);
  ps.join(pr3);

  rect_2_polygon(Point_2(4,1), Point_2(4,5), Point_2(6,5), Point_2(6,1), pr4);
  ps.join(pr4);

  // Output:
  std::list<Polygon_with_holes_2> result;
  ps.polygons_with_holes(std::back_inserter(result));

  std::copy(result.begin(), result.end(),       // export to standard output
            std::ostream_iterator<Polygon_with_holes_2>(std::cout, "\n"));
  std::cout << std::endl;
  
  return 0;
}
