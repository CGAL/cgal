//! \file examples/Arrangement_on_surface_2/rational_functions.cpp
// Constructing an arrangement of arcs of rational functions.
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
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_rational_arc_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_rational_arc_traits_2<Alg_kernel,
                                        Nt_traits>    Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Rational_arc_2;
typedef Traits_2::Rat_vector                          Rat_vector;
typedef std::list<Rational_arc_2>                     Rat_arcs_list;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main ()
{

  Traits_2 m_traits;

  // Create an arc supported by the polynomial y = 1
  Rat_vector        P1(1);
  P1[0] = 1;

  Rat_vector        P2(3);  
  P2[2] = -1; P2[1] = 0; P2[0] = 1;

  Rational_arc_2    a1 (P1, P2, Algebraic(-1), Algebraic(1));

  Rational_arc_2    a1left (P1, P2, Algebraic(-1), Algebraic(0));

  Rational_arc_2    a1right (P1, P2, Algebraic(0), Algebraic(1));

  Rational_arc_2 cv1,cv2;

  Point_2 split_point1(0,1);

  m_traits.split_2_object()(a1, split_point1, cv1, cv2);


  std::cout <<"a1left " << a1left << std::endl;
  std::cout <<"a1right " << a1right << std::endl;
  std::cout <<"cv1 " << cv1 << std::endl;
  std::cout <<"cv2 " << cv2 << std::endl;
  std::cout <<"m_traits.equal_2_object() " << m_traits.equal_2_object()(cv1,a1left) << std::endl;
  std::cout <<"m_traits.equal_2_object() " << m_traits.equal_2_object()(cv2,a1right) << std::endl;

  //***********************************************//
/*
  Rational_arc_2    a2 (P1, P2, Algebraic(-1), false);

  Point_2 split_point2(-3,-0.125);

  m_traits.split_2_object()(a2, split_point2, cv1, cv2);

  std::cout <<"cv1 " << cv1 << std::endl;
  std::cout <<"cv2 " << cv2 << std::endl;

  Rational_arc_2    a3 (P1, P2, Algebraic(1), true);

  Point_2 split_point3(3,-0.125);

  m_traits.split_2_object()(a3, split_point3, cv1, cv2);

  std::cout <<"cv1 " << cv1 << std::endl;
  std::cout <<"cv2 " << cv2 << std::endl;
*/
  return 0;
}

#endif
