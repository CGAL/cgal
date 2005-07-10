// file: examples/Arrangement_2/example17.C


//#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

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
typedef CGAL::Arr_naive_point_location<Arrangement_2> Point_location;

int main ()
{
  Arrangement_2    arr;
  Point_location   pl (arr);

  Rat_point_2      c1 = Rat_point_2 (0, 0);
  Rational         sqr_r1 = Rational (25);       // = 5^2
  Rat_circle_2     circ1 = Rat_circle_2 (c1, sqr_r1, CGAL::CLOCKWISE);
  Conic_arc_2      cv1 = Conic_arc_2 (circ1);

  Rat_point_2      c2 = Rat_point_2 (7, 7);
  Rational         sqr_r2 = Rational (25);       // = 5^2
  Rat_circle_2     circ2 = Rat_circle_2 (c2, sqr_r2, CGAL::CLOCKWISE);
  Conic_arc_2      cv2 = Conic_arc_2 (circ2);

  Rat_point_2      c3 = Rat_point_2 (4, Rational (-1,2));
  Rational         sqr_r3 = Rational (49, 4);    // = 3.5^2
  Rat_circle_2     circ3 = Rat_circle_2 (c3, sqr_r3, CGAL::CLOCKWISE);
  Conic_arc_2      cv3 = Conic_arc_2 (circ3);

  insert (arr, pl, cv1);
  insert (arr, pl, cv2);
  insert (arr, pl, cv3);
  
  // Print the arrangement vertices.
  Arrangement_2::Vertex_const_iterator  vit;
  int                                   i;

  std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
  for (i = 1, vit = arr.vertices_begin(); 
       vit != arr.vertices_end(); vit++, i++)
  {
    std::cout << '\t' << i << ": " << vit->point() << std::endl;
  }
  std::cout << std::endl;

  // Print the arrangement edges.
  Arrangement_2::Edge_const_iterator    eit;

  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  for (i = 1, eit = arr.edges_begin(); 
       eit != arr.edges_end(); eit++, i++)
  {
    std::cout << '\t' << i << ": " << eit->curve() << std::endl;
  }

  std::cout << std::endl;

  std::cout << arr.number_of_faces() << " faces." << std::endl;

  return (0);
}

