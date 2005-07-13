// file: examples/Arrangement_2/example14.C

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef Rat_kernel::Point_2                             Rat_point_2;
typedef Rat_kernel::Segment_2                           Rat_segment_2;
typedef Rat_kernel::Circle_2                            Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, 
                                 Alg_kernel,
                                 Nt_traits>             Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Conic_arc_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2>   Naive_pl;

int main ()
{
  Arrangement_2    arr;
  Naive_pl         pl (arr);

  // Insert a hyperbolic arc, supported by the hyperbola y = 1/x
  // (or: xy - 1 = 0) with the end-points (1/5, 4) and (2, 1/2).
  // Note that the arc is counterclockwise oriented.
  Point_2       ps1 (Rational(1,4), 4);
  Point_2       pt1 (2, Rational(1,2));
  Conic_arc_2   c1 (0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE, ps1, pt1);

  insert (arr, c1, pl);

  // Insert a full ellipse, which is (x/4)^2 + (y/2)^2 = 0 rotated by
  // phi=36.87 degree (such that sin(phi) = 0.6, cos(phi) = 0.8),
  // yielding: 58x^2 + 72y^2 - 48xy - 360 = 0.
  Conic_arc_2   c2 (58, 72, -48, 0, 0, -360);
  
  insert (arr, c2, pl);

  // Insert the segment (1, 1) -- (0, -3).
  Rat_point_2   ps3 (1, 1);
  Rat_point_2   pt3 (0, -3);
  Conic_arc_2   c3 (Rat_segment_2 (ps3, pt3));

  insert (arr, c3, pl);

  // Insert a circular arc supported by the circle x^2 + y^2 = 5^2,
  // with (-3, 4) and (4, 3) as its endpoints. We want the arc to be
  // clockwise oriented, so it passes through (0, 5) as well.
  Rat_point_2   ps4 (-3, 4);
  Rat_point_2   pm4 (0, 5);
  Rat_point_2   pt4 (4, 3);
  Conic_arc_2   c4 (ps4, pm4, pt4);

  CGAL_assertion (c4.is_valid());
  insert (arr, c4, pl);

  // Insert a full unit circle that is centered at (0, 4).
  Rat_circle_2  circ5 (Rat_point_2(0,4), 1);
  Conic_arc_2   c5 (circ5);
  
  insert (arr, c5, pl);

  // Insert a parabolic arc that is supported by a parabola y = -x^2
  // (or: x^2 + y = 0) and whose end-points are (-sqrt(3), -3) ~ (-1.73, -3)
  // and (sqrt(2), -2) ~ (1.41, -2). Notice that since the x-coordinates 
  // of the end-points cannot be acccurately represented, we specify them
  // as the intersections of the parabola with the lines y = -3 and y = -2.
  // Note that the arc is clockwise oriented.
  Conic_arc_2   c6 =
    Conic_arc_2 (1, 0, 0, 0, 1, 0,       // The parabola.
                 CGAL::CLOCKWISE,
                 Point_2 (-1.73, -3),    // Approximation of the source.
                 0, 0, 0, 0, 1, 3,       // The line: y = -3.
                 Point_2 (1.41, -2),     // Approximation of the target.
                 0, 0, 0, 0, 1, 2);      // The line: y = -2.

  CGAL_assertion (c6.is_valid());
  insert (arr, c6, pl);

  // Insert the right half of the circle centered at (4, 2.5) whose radius
  // is 1/2 (therefore its squared radius is 1/4).
  Rat_circle_2  circ7 (Rat_point_2(4, Rational(5,2)), Rational(1,4));
  Point_2       ps7 (4, 3);
  Point_2       pt7 (4, 2);
  Conic_arc_2   c7 (circ7, CGAL::CLOCKWISE, ps7, pt7);
  
  insert (arr, c7, pl);

  // Print out the size of the resulting arrangement.
  std::cout << "V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() << std::endl;

  return (0);
}

