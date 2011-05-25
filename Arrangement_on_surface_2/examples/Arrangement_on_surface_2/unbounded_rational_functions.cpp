//! \file examples/Arrangement_2/unbounded_rational_functions.cpp
// Constructing an arrangement of unbounded portions of rational functions.
#include <CGAL/basic.h>

#include <CGAL/Gmpz.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Gmpz                                    Integer;
typedef CGAL::Algebraic_kernel_d_1<Integer>           AK1; 
typedef CGAL::Arr_rational_function_traits_2<AK1>     Traits_2; 

typedef AK1::Polynomial_1                             Polynomial_1;
typedef AK1::Algebraic_real_1                         Alg_real_1;

typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;

int main ()
{
  // an object of Algebraic_kernel_d_1
  AK1 ak1;
  
  // Traits class object 
  Traits_2 traits(ak1); 
  
  // constructor for rational functions 
  Traits_2::Construct_curve_2 construct
    = traits.construct_curve_2_object(); 
  
  // a polynomial representing x .-)
  Polynomial_1 x = CGAL::shift(Polynomial_1(1),1);
  
  // container storing all arcs 
  std::vector<Traits_2::Curve_2>  arcs;

  
  // Create the rational functions (y = 1 / x), and (y = -1 / x).
  Polynomial_1 P1(1);
  Polynomial_1 Q1 = x;
  arcs.push_back(construct( P1, Q1));
  arcs.push_back(construct(-P1, Q1));

  // Create a bounded segments of the parabolas (y = -4*x^2 + 3) and
  // (y = 4*x^2 - 3), defined over [-sqrt(3)/2, sqrt(3)/2].
  Polynomial_1 P2 = -4*x*x+3; 
  std::vector<Alg_real_1> roots;
  ak1.solve_1_object()(P2,std::back_inserter(roots));// [-sqrt(3)/2, sqrt(3)/2]
  arcs.push_back(construct( P2, roots[0], roots[1]));
  arcs.push_back(construct(-P2, roots[0], roots[1]));

  // Create the rational function (y = 1 / 2*x) for x > 0, and the
  // rational function (y = -1 / 2*x) for x < 0.
  Polynomial_1 P3(1);
  Polynomial_1 Q3 = 2*x;
  arcs.push_back(construct( P3, Q3, Algebraic(0), true ));
  arcs.push_back(construct(-P3, Q3, Algebraic(0), false));

  // Construct the arrangement of the six arcs.
  Arrangement_2              arr;
  insert (arr, arcs.begin(), arcs.end());

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << arr.number_of_vertices()
            << " (plus " << arr.number_of_vertices_at_infinity()
            << " at infinity)"
            << ",  E = " << arr.number_of_edges() 
            << ",  F = " << arr.number_of_faces() 
            << " (" << arr.number_of_unbounded_faces() << " unbounded)"
            << std::endl << std::endl;

  return 0;
}

#endif
