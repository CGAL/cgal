
// Constructing an arrangement of polycurves.

#include <CGAL/basic.h>
#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl; 
  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <vector>
#include <list>
#include "arr_print.h"

#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
//#include <CGAL/Arr_algebraic_segment_traits_2.h>


typedef CGAL::Quotient<CGAL::MP_Float>                                    Number_type;
typedef CGAL::Cartesian<Number_type>                                      Kernel;

//////////////
//line segment traits
//////////////
typedef CGAL::Arr_segment_traits_2<Kernel>                                Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>                     Polycurve_segment_traits_2;
typedef CGAL::Arrangement_2<Polycurve_segment_traits_2>                   Segment_arrangment_2;
typedef Polycurve_segment_traits_2::Point_2                               Segment_point_2;
typedef Polycurve_segment_traits_2::Curve_2                               Segment_section_2;

///////////////
//circle segment traits
//////////////
typedef CGAL::Arr_circle_segment_traits_2<Kernel>                         Arc_traits_2;
typedef CGAL::Arr_polyline_traits_2<Arc_traits_2>                         Polycurve_arc_traits_2;
typedef Arc_traits_2::CoordNT                                             CoordNT;
typedef Arc_traits_2::Point_2                                             Arc_point_2;
typedef Arc_traits_2::Curve_2                                             Arc_section_2;
typedef CGAL::Arrangement_2<Polycurve_arc_traits_2>                       Arc_arrangment_2;

//to be removed later
typedef CGAL::Arrangement_2<Arc_traits_2>                 Arrangement_2_test;

////////////////////
//conic traits
////////////////////
typedef CGAL::CORE_algebraic_number_traits                                Nt_traits;
typedef Nt_traits::Rational                                               Rational;
typedef Nt_traits::Algebraic                                              Algebraic;
typedef CGAL::Cartesian<Rational>                                         Rat_kernel;
typedef Rat_kernel::Point_2                                               Rat_point_2;
typedef Rat_kernel::Segment_2                                             Rat_segment_2;
typedef Rat_kernel::Circle_2                                              Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                                        Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>       Conic_traits_2;
typedef Conic_traits_2::Point_2                                           Conic_point_2;
typedef Conic_traits_2::Curve_2                                           Conic_arc_2;
typedef CGAL::Arr_polyline_traits_2<Conic_traits_2>                       Polycurve_conic_traits_2;
typedef CGAL::Arrangement_2<Polycurve_conic_traits_2>                     Conic_arrangment_2; 

///////////////////
//Bezier Curves
///////////////////
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>    Bezier_traits_2;
typedef CGAL::Arr_polyline_traits_2<Bezier_traits_2>                          Polycurve_bezier_traits_2;
typedef CGAL::Arrangement_2<Polycurve_bezier_traits_2>                        Bezier_arrangment_2;
typedef Bezier_traits_2::Point_2                                              Bezier_point_2;
typedef Bezier_traits_2::Curve_2                                              Bezier_section_2;  

//////////////////
//Algebric Curves
//////////////////
// #if CGAL_USE_GMP && CGAL_USE_MPFI
// #include <CGAL/Gmpz.h>
// typedef CGAL::Gmpz Integer;
// #elif CGAL_USE_CORE
// #include <CGAL/CORE_BigInt.h>
// typedef CORE::BigInt Integer;
// #else
// #include <CGAL/leda_integer.h>
// typedef LEDA::integer Integer;
// #endif

// typedef CGAL::Arr_algebraic_segment_traits_2<Integer>                     Algebric_traits_2;
// typedef Algebric_traits_2::Curve_2                                        Algebric_section_2;
// typedef Algebric_traits_2::Polynomial_2                                   Algebric_polynomial_2;
// typedef CGAL::Arr_polyline_traits_2<Algebric_traits_2>                    Polycurve_algebric_traits_2;
// typedef CGAL::Arrangement_2<Polycurve_algebric_traits_2>                  Algebric_arrangment_2;




void Make_linear_polycurve (Segment_arrangment_2* Curve_arr)
{

  Segment_point_2               points1[5];
  points1[0] = Segment_point_2 (0, 0);
  points1[1] = Segment_point_2 (2, 4);
  points1[2] = Segment_point_2 (3, 0);
  points1[3] = Segment_point_2 (4, 4);
  points1[4] = Segment_point_2 (6, 0);
  Segment_section_2            pi1 (&points1[0], &points1[5]);

  std::list<Segment_point_2>    points2;
  points2.push_back (Segment_point_2 (1, 3));
  points2.push_back (Segment_point_2 (0, 2));
  points2.push_back (Segment_point_2 (1, 0));
  points2.push_back (Segment_point_2 (2, 1));
  points2.push_back (Segment_point_2 (3, 0));
  points2.push_back (Segment_point_2 (4, 1));
  points2.push_back (Segment_point_2 (5, 0));
  points2.push_back (Segment_point_2 (6, 2));
  points2.push_back (Segment_point_2 (5, 3));
  points2.push_back (Segment_point_2 (4, 2));
  Segment_section_2            pi2 (points2.begin(), points2.end());

  std::vector<Segment_point_2>  points3 (4);
  points3[0] = Segment_point_2 (0, 2);
  points3[1] = Segment_point_2 (1, 2);
  points3[2] = Segment_point_2 (3, 6);
  points3[3] = Segment_point_2 (5, 2);
  Segment_section_2            pi3 (points3.begin(), points3.end());
  
  insert (*Curve_arr, pi1);
  insert (*Curve_arr, pi2);
  insert (*Curve_arr, pi3);
}

void Make_conic_polycurve(Conic_arrangment_2* Conic_arr)
{

  // Insert a hyperbolic arc, supported by the hyperbola y = 1/x
  // (or: xy - 1 = 0) with the endpoints (1/5, 4) and (2, 1/2).
  // Note that the arc is counterclockwise oriented.
  Conic_point_2       ps1 (Rational(1,4), 4);
  Conic_point_2       pt1 (2, Rational(1,2));
  Conic_arc_2         c1 (0, 0, 1, 0, 0, -1, CGAL::COUNTERCLOCKWISE, ps1, pt1);

  // Insert a parabolic arc that is supported by a parabola y = -x^2
  // (or: x^2 + y = 0) and whose endpoints are (-sqrt(3), -3) ~ (-1.73, -3)
  // and (sqrt(2), -2) ~ (1.41, -2). Notice that since the x-coordinates 
  // of the endpoints cannot be acccurately represented, we specify them
  // as the intersections of the parabola with the lines y = -3 and y = -2.
  // Note that the arc is clockwise oriented.
  Conic_arc_2   c2 =
    Conic_arc_2 (1, 0, 0, 0, 1, 0,       // The parabola.
                 CGAL::CLOCKWISE,
                 Conic_point_2 (-1.73, -3),    // Approximation of the source.
                 0, 0, 0, 0, 1, 3,       // The line: y = -3.
                 Conic_point_2 (1.41, -2),     // Approximation of the target.
                 0, 0, 0, 0, 1, 2);      // The line: y = -2.
  CGAL_assertion (c2.is_valid());
  
  // Insert the segment (1, 1) -- (0, -3).
  Rat_point_2   ps3 (1, 1);
  Rat_point_2   pt3 (0, -3);
  Conic_arc_2   c3 (Rat_segment_2 (ps3, pt3));
  
   insert (*Conic_arr, c1);
   insert (*Conic_arr, c2);
   insert (*Conic_arr, c3);
}

void Make_arc_polycurve( Arc_arrangment_2* Arc_arr)
{
  std::list<Arc_section_2>  curves;

  // Create a circular arc that correspond to the upper half of the
  // circle centered at (1,1) with squared radius 3. We create the
  // circle with clockwise orientation, so the arc is directed from
  // (1 - sqrt(3), 1) to (1 + sqrt(3), 1).
  Kernel::Point_2   c1 = Kernel::Point_2 (1, 1);
  Kernel::Circle_2  circ1 = Kernel::Circle_2 (c1, 3, CGAL::CLOCKWISE);
  CoordNT           one_minus_sqrt_3 = CoordNT (1, -1, 3);
  CoordNT           one_plus_sqrt_3 = CoordNT (1, 1, 3);
  Arc_point_2       s1 = Arc_point_2 (one_minus_sqrt_3, CoordNT (1));
  Arc_point_2       t1 = Arc_point_2 (one_plus_sqrt_3, CoordNT (1));
  curves.push_back (Arc_section_2 (circ1, s1, t1));

  //insert ( *Arc_arr, curves.begin(), curves.end());

  //to be removed later
  //Arrangement_2_test arr;
  //insert(arr, curves.begin(), curves.end() );
}

void Make_bezier_polycurve ( Bezier_arrangment_2* bezier_arr) 
{
 
  const char* filename = "Bezier.dat";

  // Open the input file.
  std::ifstream  in_file (filename);

  if (! in_file.is_open()) 
  {
    std::cerr << "Failed to open " << filename << std::endl;
    return ;
  }

  // Read the curves from the input file.
  unsigned int                  n_curves;
  std::list<Bezier_section_2>   curves;
  Bezier_section_2              B;
  unsigned int                  k;
  
  in_file >> n_curves;
  
  for (k = 0; k < n_curves; k++) 
  {
    // Read the current curve (specified by its control points).
    in_file >> B;
    curves.push_back (B);
  }

  //insert ( *bezier_arr, curves.begin(), curves.end());

  // Construct the arrangement for testing.
  CGAL::Arrangement_2<Bezier_traits_2>    arr;
  insert( arr, curves.begin(), curves.end() );
}

void Make_algebric_polycurve()
{
  //TODO: 
}

int main ()
{
 
  Segment_arrangment_2        Linear_segment_arr;
  Arc_arrangment_2            Arc_arr;
  Conic_arrangment_2          Conic_arr;
  Bezier_arrangment_2         Bezier_arr;
  //Algebric_arrangment_2       Algebric_arr;
  
  Make_linear_polycurve ( & Linear_segment_arr );
  Make_conic_polycurve  ( & Conic_arr          );
  Make_arc_polycurve    ( & Arc_arr            );
  Make_bezier_polycurve ( & Bezier_arr         );
  
  //print_arrangement (Linear_segment_arr);
  //print_arrangement (Conic_arr);

  
  return 0;
}
#endif
