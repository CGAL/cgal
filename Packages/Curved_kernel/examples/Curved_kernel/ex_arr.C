// inspired from file: examples/Arrangement_2/example4.C


//#include "short_names.h"

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Algebraic_kernel_2_2.h>

#include <CGAL/NT_extensions_Root_of/CGAL_Quotient.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Gmpq.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Lazy_exact_nt.h>

#include <CGAL/Circular_kernel.h>
#include <CGAL/Circular_arc_traits.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>

#include <CGAL/Circular_arc_traits_checker.h>

// typedef CGAL::MP_Float                                      NT;
typedef CGAL::Cartesian<CGAL::MP_Float>                     Linear_k;

typedef CGAL::Algebraic_kernel_2_2<Linear_k::RT>            Algebraic_k;
typedef CGAL::Curved_kernel<Linear_k,Algebraic_k>           Curved_k;
typedef CGAL::Circular_arc_traits<Curved_k>                 T0;

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef Rat_kernel::Point_2                           Rat_point_1;
typedef Rat_kernel::Circle_2                          Rat_circle_1;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, 
				 Alg_kernel,
				 Nt_traits>           T1;

typedef CGAL::Circular_arc_traits_checker<T0, T1>     T2;

// to enter input data
typedef Linear_k::Point_2                           Rat_point_0;
typedef Linear_k::Circle_2                          Rat_circle_0;
typedef Linear_k::RT                                NT;
//

typedef T0::Point_2                                   Point_0;
typedef T0::Curve_2                                   Curve_0;
typedef T1::Point_2                                   Point_1;
typedef T1::Curve_2                                   Curve_1;

typedef T2::Point_2                                   Point_2;
typedef T2::Curve_2                                   Curve_2;
typedef CGAL::Arrangement_2<T2>                       Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2> Point_location;

int main ()
{
  Arrangement_2    arr;
  Point_location   pl (arr);

  Rat_point_0      c1_0 = Rat_point_0 (0, 0);
  NT               sqr_r1_0 = NT (25);       // = 5^2
  Rat_circle_0     circ1_0 = Rat_circle_0 (c1_0, sqr_r1_0/*, CGAL::CLOCKWISE*/);
  Curve_0          cv1_0 = Curve_0 (circ1_0);

  Rat_point_0      c2_0 = Rat_point_0 (7, 7);
  NT               sqr_r2_0 = NT (25);       // = 5^2
  Rat_circle_0     circ2_0 = Rat_circle_0 (c2_0, sqr_r2_0/*, CGAL::CLOCKWISE*/);
  Curve_0          cv2_0 = Curve_0 (circ2_0);

  Rat_point_0      c3_0 = Rat_point_0 (4, NT (-0.5)); // was (-1,2)
  NT               sqr_r3_0 = NT (12.25);    // = 3.5^2 was (49, 4)
  Rat_circle_0     circ3_0 = Rat_circle_0 (c3_0, sqr_r3_0/*, CGAL::CLOCKWISE*/);
  Curve_0          cv3_0 = Curve_0 (circ3_0);


  Rat_point_1      c1_1 = Rat_point_1 (0, 0);
  Rational         sqr_r1_1 = Rational (25);       // = 5^2
  Rat_circle_1     circ1_1 = Rat_circle_1 (c1_1, sqr_r1_1, CGAL::CLOCKWISE);
  Curve_1          cv1_1 = Curve_1 (circ1_1);

  Rat_point_1      c2_1 = Rat_point_1 (7, 7);
  Rational         sqr_r2_1 = Rational (25);       // = 5^2
  Rat_circle_1     circ2_1 = Rat_circle_1 (c2_1, sqr_r2_1, CGAL::CLOCKWISE);
  Curve_1          cv2_1 = Curve_1 (circ2_1);

  Rat_point_1      c3_1 = Rat_point_1 (4, Rational (-1,2));
  Rational         sqr_r3_1 = Rational (49, 4);    // = 3.5^2
  Rat_circle_1     circ3_1 = Rat_circle_1 (c3_1, sqr_r3_1, CGAL::CLOCKWISE);
  Curve_1          cv3_1 = Curve_1 (circ3_1);

  std::cout << std::endl << "insertion cercle " 
	    << "_1 " << cv1_0 << std::endl
	    << "_2 " << cv1_1 << std::endl << std::endl;
  insert (arr, pl, Curve_2(cv1_0, cv1_1));

  std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  std::cout << arr.number_of_faces() << " faces." << std::endl;

  std::cout << std::endl << "insertion cercle " 
	    << "_1 " << cv2_0 << std::endl << std::endl
	    << "_2 " << cv2_1 << std::endl;
  insert (arr, pl, Curve_2(cv2_0, cv2_1));

  std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  std::cout << arr.number_of_faces() << " faces." << std::endl;

  std::cout << std::endl << "insertion cercle " 
	    << "_1 " << cv3_0 << std::endl
	    << "_2 " << cv3_1 << std::endl << std::endl;
  insert (arr, pl, Curve_2(cv3_0, cv3_1));
  
  // Print the arrangement vertices.
  Arrangement_2::Vertex_const_iterator  vit;
  Arrangement_2::Vertex_const_handle    vh;
  int                                   i;

  std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
  for (i = 1, vit = arr.vertices_begin(); 
       vit != arr.vertices_end(); vit++, i++)
  {
    vh = *vit;
    std::cout << '\t' << i << ": pair < " << vh.point().first << " , " 
	      << vh.point().second << " > " << std::endl;
  }
  std::cout << std::endl;

  // Print the arrangement edges.
  Arrangement_2::Edge_const_iterator    eit;
  Arrangement_2::Halfedge_const_handle  hh;

  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  for (i = 1, eit = arr.edges_begin(); eit != arr.edges_end(); eit++, i++)
  {
    hh = *eit;
    std::cout << '\t' << i << ": pair < " << hh.curve().first << " , "
	      << hh.curve().second << " > " << std::endl;
  }
  std::cout << std::endl;

  std::cout << arr.number_of_faces() << " faces." << std::endl;

  return (0);
}

