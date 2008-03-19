// Test for CGAL::constant<T,i>()
// Sylvain Pion, 2006.

#include <CGAL/constant.h>
#include <CGAL/MP_Float.h>
#include <cassert>

// Just to have a look at the assembly code if the compiler
// can propagate constants through it.
// g++ 3.4 cannot, but g++ 4.0 and 4.1 can.
template < typename T >
T cprop()
{
  return CGAL::constant<double, 2>()
       + CGAL::constant<double, 3>()
       + CGAL::constant<double, 4>();
}

// Hard to check that the function returns a const-ref,
// however we can check that the return address is always the same.
template < typename T >
void
check_address()
{
  const T * ptr1 = & CGAL::constant<T, 1>();
  const T * ptr2 = & CGAL::constant<T, 1>();
  assert( ptr1 == ptr2 );
}

int main()
{
  const double & d = CGAL::constant<double, 2>();
  const CGAL::MP_Float & m = CGAL::constant<CGAL::MP_Float, -3>();

  assert( d == 2 );
  assert( m == -3 );

  check_address<double>();
  check_address<CGAL::MP_Float>();

  cprop<double>();
  cprop<CGAL::MP_Float>(); // dream...

  return 0;
}
