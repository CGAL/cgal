//! \file examples/Arrangement_2/example13.C
// Constructing an arrangement of circles using the conic-arc traits.

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef Rat_kernel::Point_2                           Rat_point_2;
typedef Rat_kernel::Circle_2                          Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, 
                                 Alg_kernel,
                                 Nt_traits>           Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Conic_arc_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main ()
{
  // Create a circle centered at the origin with radius 5.
  Rat_point_2      c1 = Rat_point_2 (0, 0);
  Rational         sqr_r1 = Rational (25);       // = 5^2
  Rat_circle_2     circ1 = Rat_circle_2 (c1, sqr_r1, CGAL::CLOCKWISE);
  Conic_arc_2      cv1 = Conic_arc_2 (circ1);

  // Create a circle centered at (7,7) with radius 5.
  Rat_point_2      c2 = Rat_point_2 (7, 7);
  Rational         sqr_r2 = Rational (25);       // = 5^2
  Rat_circle_2     circ2 = Rat_circle_2 (c2, sqr_r2, CGAL::CLOCKWISE);
  Conic_arc_2      cv2 = Conic_arc_2 (circ2);

  // Create a circle centered at (4,-0.5) with radius 3.5 (= 7/2).
  Rat_point_2      c3 = Rat_point_2 (4, Rational (-1,2));
  Rational         sqr_r3 = Rational (49, 4);    // = 3.5^2
  Rat_circle_2     circ3 = Rat_circle_2 (c3, sqr_r3, CGAL::CLOCKWISE);
  Conic_arc_2      cv3 = Conic_arc_2 (circ3);

  // Construct the arrangement of the three circles.
  Arrangement_2    arr;

  insert_curve (arr, cv1);
  insert_curve (arr, cv2);
  insert_curve (arr, cv3);
  
  // Locate the vertex with maximal degree.
  Arrangement_2::Vertex_const_iterator  vit;
  Arrangement_2::Vertex_const_handle    v_max;
  unsigned int                          max_degree = 0;

  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    if (vit->degree() > max_degree)
    {
      v_max = vit;
      max_degree = vit->degree();
    }
  }

  std::cout << "v_max = (" << v_max->point() << ") "
            << "with degree " << max_degree << "." << std::endl;

  return (0);
}

