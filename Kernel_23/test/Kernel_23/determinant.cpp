#define CGAL_CHECK_DETERMINANT

#include <CGAL/Simple_cartesian.h>

int main()
{
  assert(CGAL::determinant<int>(4, 5, 1, 4, 6, 3, 1,
                                4, 3, 6, 4, 2, 7, 3,
                                6, 3, 3, 6, 2, 4, 5,
                                1, 4, 3, 5, 5, 6 ,1,
                                1, 3, 2, 7, 9, 6, 1,
                                7, 6, 5, 4, 6, 2, 2,
                                2, 3, 5, 7, 4, 3, 3) == 763);
  return 0;
}
