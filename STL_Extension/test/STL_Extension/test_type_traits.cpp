// ============================================================================
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// $URL: $
// $Id: $
// author(s)     : Andreas Meyer <ameyer@mpi-inf.mpg.de>


#include <cassert>
#include <CGAL/type_traits.h>

struct A {};
struct B : public A {};
typedef A C;

int main() {
  static_assert( ( ::CGAL::is_same_or_derived< A,A >::value == 1 ) );
  static_assert( ( ::CGAL::is_same_or_derived< A,B >::value == 1 ) );
  static_assert( ( ::CGAL::is_same_or_derived< B,A >::value == 0 ) );
  static_assert( ( ::CGAL::is_same_or_derived< B,B >::value == 1 ) );
  static_assert( ( ::CGAL::is_same_or_derived< A,C >::value == 1 ) );
  static_assert( ( ::CGAL::is_same_or_derived< B,C >::value == 0 ) );
  static_assert( ( ::CGAL::is_same_or_derived< C,C >::value == 1 ) );
  static_assert( ( ::CGAL::is_same_or_derived< C,A >::value == 1 ) );
  static_assert( ( ::CGAL::is_same_or_derived< C,B >::value == 1 ) );

  static_assert(   CGAL::is_convertible_without_narrowing_v<int, int> );
  static_assert( ! CGAL::is_convertible_without_narrowing_v<int, signed char> );
  static_assert(   CGAL::is_convertible_without_narrowing_v<signed char, int> );
  static_assert( ! CGAL::is_convertible_without_narrowing_v<int, double> );
  static_assert(   CGAL::is_convertible_without_narrowing_v<float, double> );
  static_assert( ! CGAL::is_convertible_without_narrowing_v<double, float> );
  static_assert( ! CGAL::is_convertible_without_narrowing_v<A, int> );
  return 0;
}
