#include <CGAL/known_bit_size_integers.h>
#include <CGAL/Testsuite/assert.h>
#include <iostream>

int main()
{
  std::cout << "Verifying the sizes of CGAL::[U]Integer{8,16,32,64}"
            << std::endl;

  CGAL_test_assert(sizeof(CGAL::Integer8)  == 1);
  CGAL_test_assert(sizeof(CGAL::Integer16) == 2);
  CGAL_test_assert(sizeof(CGAL::Integer32) == 4);
#ifdef CGAL_HAS_INTEGER64
  CGAL_test_assert(sizeof(CGAL::Integer64) == 8);
#endif

  CGAL_test_assert(sizeof(CGAL::UInteger8)  == 1);
  CGAL_test_assert(sizeof(CGAL::UInteger16) == 2);
  CGAL_test_assert(sizeof(CGAL::UInteger32) == 4);
#ifdef CGAL_HAS_INTEGER64
  CGAL_test_assert(sizeof(CGAL::UInteger64) == 8);
#endif

  return 0;
}
