#include <CGAL/known_bit_size_integers.h>
#include <cassert>
#include <iostream>

int main()
{
  std::cout << "Verifying the sizes of CGAL::[U]Integer{8,16,32,64}"
            << std::endl;

  assert(sizeof(CGAL::Integer8)  == 1);
  assert(sizeof(CGAL::Integer16) == 2);
  assert(sizeof(CGAL::Integer32) == 4);
#ifdef CGAL_HAS_INTEGER64
  assert(sizeof(CGAL::Integer64) == 8);
#endif

  assert(sizeof(CGAL::UInteger8)  == 1);
  assert(sizeof(CGAL::UInteger16) == 2);
  assert(sizeof(CGAL::UInteger32) == 4);
#ifdef CGAL_HAS_INTEGER64
  assert(sizeof(CGAL::UInteger64) == 8);
#endif

  return 0;
}
