#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_cached_traits_2.h>
#include <CGAL/Quotient.h>

#include "include/Segment_traits_test.h"

typedef CGAL::Quotient<int>                             NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_cached_traits_2<Kernel>       Traits;
typedef Traits::Segment_cached_2                        Segment_cached_2;

std::ostream & operator<<(std::ostream & os, const Segment_cached_2 & seg)
{
  os << static_cast<Kernel::Segment_2>(seg);
  return (os);
}

int main(int argc, char * argv[])
{
  Segment_traits_test< Traits, NT >  test_obj( argc, argv );
  return (test_obj.start()) ? 0 /* SUCCESS */ : 1; /* FAILURE */
}
