/* Test the insertion of a single curve that has 2 vertical asymptotes.
 *
 *      2      5
 * o----o------o----o
 * |    |      |    |
 * |    \      /    |
 * |     ------     |
 * |     ------     |
 * |    /      \    |
 * |    |      |    |
 * o----------------o
 *
 * When the maximal end of the curve is processed, a vertex that corresponds
 * to the minimal end has already been inserted into the arrangement, but the
 * pair of halfedges that correspond to the curve hasn't been created yet.
 * The vertex has degree 2, and the arrangement is in an intermediate state.
 * Thus, a special treatment is required.
 */

#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
//  bool   UNTESTED_TRAITS_AS_CORE_IS_NOT_ISTALLED;
  std::cout << std::endl
            << "NOTE: WARNING: Core is not installed, "
            << "skipping the test of the conic traits ..."
            << std::endl;
  return 0;
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_rational_arc_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::CORE_algebraic_number_traits                      Nt_traits;
typedef Nt_traits::Rational                                     Rational;
typedef Nt_traits::Algebraic                                    Algebraic;
typedef CGAL::Cartesian<Algebraic>                              Alg_kernel;
typedef CGAL::Arr_rational_arc_traits_2<Alg_kernel, Nt_traits>  Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef Traits_2::Curve_2                                       Rational_arc_2;
typedef Traits_2::Rat_vector                                    Rat_vector;
typedef std::list<Rational_arc_2>                               Rat_arcs_list;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;

int main ()
{
  std::list<Rational_arc_2>  arcs;

  // Create the rational function y = 1 / ((x - 2)(x - 5)) = 1 / (x^2 - 7x + 10)
  // Create the rational function y = 1 / ((2 - x)(x - 5)) = 1 / (-x^2 + 7x - 10)
  Rat_vector        P1(1);
  P1[0] = 1;

  Rat_vector        Q1(3);
  Rat_vector        Q2(3);
  Q1[2] = 1; Q1[1] = -7; Q1[0] = 10;
  Q2[2] = -1; Q2[1] = 7; Q2[0] = -10;

  Rational_arc_2 c1(P1, Q1, Algebraic(2), Algebraic(5));
  Rational_arc_2 c2(P1, Q2, Algebraic(2), Algebraic(5));

  // Construct the arrangement of the six arcs.
  Arrangement_2 arr;
  insert (arr, c1);
  insert (arr, c2);
  
  if (!arr.is_valid()) {
    std::cerr << "The arrangement is not valid!" << std::endl;
    return -1;
  }
  
  const char * names[5] = {
    "number of vertices",
    "number of vertices at infinity",
    "number of edges",
    "number of faces",
    "number of unbounded faces"
  };
  
  Arrangement_2::Size expected_sizes[] = {0, 4, 2, 3, 3};
  Arrangement_2::Size sizes[5];
  sizes[0] = arr.number_of_vertices();
  sizes[1] = arr.number_of_vertices_at_infinity();
  sizes[2] = arr.number_of_edges();
  sizes[3] = arr.number_of_faces();
  sizes[4] = arr.number_of_unbounded_faces();

  unsigned int i;
  int result = 0;
  for (i = 0; i < 5; ++ i) {
    if (expected_sizes[i] != sizes[i]) {
      std::cerr << names[i] << ": " << sizes[i] << ", expected: "
                << expected_sizes[i] << std::endl;
      result = -1;
    } else {
      std::cout << names[i] << ": " << sizes[i] << std::endl;
    }
  }

  return result;
}

#endif
