#define TEST_CARTESIAN 1
#define KERNEL CGAL::Cartesian< FT >
#define LAZY_KERNEL CGAL::Cartesian< LFT >

#include "test_periodic_3_triangulation_traits_3.h"

int main(int, char**)
{
  test_periodic_3_triangulation_traits_3();
}
