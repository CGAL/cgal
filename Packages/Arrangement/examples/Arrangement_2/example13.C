// examples/Arrangement_2/example3.C   
// ---------------------------------

#include "short_names.h"

#ifndef CGAL_USE_LEDA
// To enable compilation without leda:
int main ()
{
  return (0);
}

#else

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/leda_real.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_conic_traits_2.h> 
#include <CGAL/Arrangement_2.h>

typedef leda_real                                       NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_conic_traits_2<Kernel>                Traits;
typedef Traits::Point_2                                 Point_2;
typedef Traits::Segment_2                               Segment_2;
typedef Traits::Circle_2                                Circle_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef CGAL::Arr_base_node<Curve_2, X_monotone_curve_2> Base_node;
typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits,Base_node>      Arr_2;

int main()
{
  Arr_2 arr;  
             
  // Insert a hyperbolic arc, supported by the hyperbola y = 1/x
  // (or: xy - 1 = 0) with the end-points (0.25, 4) and (2, 0.5).
  Point_2   ps1 (0.25, 4);
  Point_2   pt1 (2, 0.5);
  Curve_2   c1 (0, 0, 1, 0, 0, -1, ps1, pt1);

  arr.insert(c1);

  // Insert a full ellipse, which is (x/4)^2 + (y/2)^2 = 0 rotated by
  // phi=36.87 degree (such that sin(phi) = 0.6, cos(phi) = 0.8),
  // yielding: 58x^2 + 72y^2 - 48xy - 360 = 0.
  Curve_2   c2 (58, 72, -48, 0, 0, -360);
  
  arr.insert(c2);

  // Insert the segment (1, 1) -- (0, -3).
  Point_2   ps3 (1, 1);
  Point_2   pt3 (0, -3);
  Curve_2   c3 (Segment_2 (ps3, pt3));

  arr.insert(c3);

  // Insert a circular arc supported by the circle x^2 + y^2 = 5^2,
  // with (-3, 4) and (4, 3) as its endpoints.
  Point_2   ps4 (-3, 4);
  Point_2   pt4 (4, 3);
  Circle_2  circ4 (Point_2(0,0), 5*5, CGAL::CLOCKWISE);
  Curve_2   c4 (circ4, ps4, pt4);

  arr.insert(c4);

  // Insert a full unit circle that is centered at (0, 4).
  Circle_2  circ5 (Point_2(0,4), 1*1, CGAL::COUNTERCLOCKWISE);
  Curve_2   c5 (circ5);
  
  arr.insert(c5);

  // Insert a parabolic arc that is supported by a parabola y = -x^2
  // (or: x^2 + y = 0) and whose end-points are (-sqrt(3), -3) ~ (-1.73, -3)
  // and (sqrt(2), -2) ~ (1.41, -2). Notice that since the x-coordinates 
  // of the end-points cannot be acccurately represented, we specify them
  // as the intersections of the parabola with the lines y = -3 and y = -2.
  Curve_2   c6 (1, 0, 0, 0, 1, 0,       // The parabola.
		Point_2 (-1.73, -3),    // Approximation of the source.
		0, 0, 0, 0, 1, 3,       // The line: y = -3.
		Point_2 (1.41, -2),     // Approximation of the target.
		0, 0, 0, 0, 1, 2);      // The line: y = -2.

  arr.insert(c6);

  // Print out the number of vertices, edges and faces in the arrangement.
  std::cout << "Number of vertices: " 
	    << arr.number_of_vertices() << std::endl;
  
  std::cout << "Number of edges: " 
	    << arr.number_of_halfedges()/2 << std::endl;
  std::cout << "Number of faces: " 
	    << arr.number_of_faces() << std::endl;

  return 0;
}

#endif
