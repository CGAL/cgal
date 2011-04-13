#include <CGAL/_test_sign.h>
#ifdef USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_real.h>
#include <CGAL/leda_rational.h>
#endif // USE_LEDA

using namespace CGAL;

int
main()
{
#ifdef USE_LEDA
  assert( _test_sign( leda_integer(1)));
  assert( _test_sign( leda_rational(1)));
  assert( _test_sign( leda_real(1)));
#endif // USE_LEDA

  return 0;
}

