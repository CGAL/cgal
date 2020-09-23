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

#include <CGAL/config.h>

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

#include <CGAL/CORE_BigInt.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CORE::BigInt                                   Number_type;
typedef CGAL::Algebraic_kernel_d_1<Number_type>           AK1;
typedef CGAL::Arr_rational_function_traits_2<AK1>  Traits_2;

typedef Traits_2::Curve_2                          Curve_2;
typedef Traits_2::Polynomial_1                     Polynomial_1;
typedef Traits_2::Algebraic_real_1                 Alg_real_1;

typedef CGAL::Arrangement_2<Traits_2>              Arrangement_2;

int main()
{
  // Traits class object
  AK1 ak1;
  Traits_2 traits(&ak1);

  // constructor for rational functions
  Traits_2::Construct_curve_2 construct = traits.construct_curve_2_object();

  // a polynomial representing x .-)
  Polynomial_1 x = CGAL::shift(Polynomial_1(1),1);

  // Create the rational function y = 1 / ((x - 2)(x - 5)) = 1 / (x^2 - 7x + 10)
  Polynomial_1 P1(1);
  Polynomial_1 P2(-1);
  Polynomial_1 Q = x*x - 7*x + 10;

  Curve_2 c1 = construct(P1, Q, Alg_real_1(2), Alg_real_1(5));
  Curve_2 c2 = construct(P2, Q, Alg_real_1(2), Alg_real_1(5));

  // Output
  const char * names[5] = {
    "number of vertices",
    "number of vertices at infinity",
    "number of edges",
    "number of faces",
    "number of unbounded faces"
  };

  Arrangement_2::Size expected_sizes[] = {0, 2, 1, 2, 2};

  // Construct the 1st arrangement.
  Arrangement_2 arr1(&traits);
  arr1.insert_in_face_interior(c1, arr1.reference_face());

  if (!arr1.is_valid()) {
    std::cerr << "The first arrangement is not valid!" << std::endl;
    return -1;
  }

  Arrangement_2::Size sizes[5];
  sizes[0] = arr1.number_of_vertices();
  sizes[1] = arr1.number_of_vertices_at_infinity();
  sizes[2] = arr1.number_of_edges();
  sizes[3] = arr1.number_of_faces();
  sizes[4] = arr1.number_of_unbounded_faces();

  std::cout << "First arrangement:" << std::endl;
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

  // If a failure has already occurred, abort.
  if (result < 0) return result;

  std::cout << std::endl;

  // Construct the 2nd arrangement.
  Arrangement_2 arr2(&traits);
  arr2.insert_in_face_interior(c2, arr2.reference_face());

  if (!arr2.is_valid()) {
    std::cerr << "The second arrangement is not valid!" << std::endl;
    return -1;
  }

  sizes[0] = arr2.number_of_vertices();
  sizes[1] = arr2.number_of_vertices_at_infinity();
  sizes[2] = arr2.number_of_edges();
  sizes[3] = arr2.number_of_faces();
  sizes[4] = arr2.number_of_unbounded_faces();

  std::cout << "Second arrangement:" << std::endl;
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
