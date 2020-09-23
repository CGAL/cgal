#ifdef CGAL_TEST_SUITE
#  undef CGAL_TEST_SUITE
#endif

#include <CGAL/Random.h>

int main()
{
  int  u = CGAL::get_default_random().get_int( 0, 1000);
  if(u >= 0 && u < 1000)
    return 0;
  else
    return 1;
}
