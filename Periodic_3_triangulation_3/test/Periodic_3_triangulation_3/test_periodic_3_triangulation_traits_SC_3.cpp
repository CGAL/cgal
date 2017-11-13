#define TEST_CARTESIAN 1
#define KERNEL CGAL::Simple_cartesian< FT >
#define LAZY_KERNEL CGAL::Simple_cartesian< LFT >

#include "test_periodic_3_triangulation_traits_3.h"

int main(int, char**)
{
  test_periodic_3_triangulation_traits_3();
}
