// This tests the rounding mode.

#include <CGAL/Interval_arithmetic.h>

int main()
{
   bool flag = true;

   flag = flag && (CGAL_FPU_get_rounding_mode() == CGAL_FPU_NEAREST);
   cout << "default: " << (int) flag << endl;

   CGAL_FPU_set_rounding_to_zero();
   flag = flag && (CGAL_FPU_get_rounding_mode() == CGAL_FPU_ZERO);
   cout << "zero   : " << (int) flag << endl;

   CGAL_FPU_set_rounding_to_infinity();
   flag = flag && (CGAL_FPU_get_rounding_mode() == CGAL_FPU_PLUS_INFINITY);
   cout << "+inf   : " << (int) flag << endl;

   CGAL_FPU_set_rounding_to_minus_infinity();
   flag = flag && (CGAL_FPU_get_rounding_mode() == CGAL_FPU_MINUS_INFINITY);
   cout << "-inf   : " << (int) flag << endl;

   CGAL_FPU_set_rounding_to_nearest();
   flag = flag && (CGAL_FPU_get_rounding_mode() == CGAL_FPU_NEAREST);
   cout << "near   : " << (int) flag << endl;

   return (int) !flag;
}
