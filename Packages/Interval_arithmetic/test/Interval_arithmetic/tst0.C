// This tests the rounding mode.

#include <CGAL/Interval_arithmetic.h>

int main()
{
   bool flag = true;

   flag = flag && (CGAL::FPU_get_rounding_mode() == CGAL::FPU_NEAREST);
   cout << "default: " << (int) flag << endl;

   // CGAL::FPU_set_rounding_to_zero();
   CGAL::FPU_set_rounding_mode(CGAL::FPU_ZERO);
   flag = flag && (CGAL::FPU_get_rounding_mode() == CGAL::FPU_ZERO);
   cout << "zero   : " << (int) flag << endl;

   // CGAL::FPU_set_rounding_to_infinity();
   CGAL::FPU_set_rounding_mode(CGAL::FPU_PLUS_INFINITY);
   flag = flag && (CGAL::FPU_get_rounding_mode() == CGAL::FPU_PLUS_INFINITY);
   cout << "+inf   : " << (int) flag << endl;

   // CGAL::FPU_set_rounding_to_minus_infinity();
   CGAL::FPU_set_rounding_mode(CGAL::FPU_MINUS_INFINITY);
   flag = flag && (CGAL::FPU_get_rounding_mode() == CGAL::FPU_MINUS_INFINITY);
   cout << "-inf   : " << (int) flag << endl;

   // CGAL::FPU_set_rounding_to_nearest();
   CGAL::FPU_set_rounding_mode(CGAL::FPU_NEAREST);
   flag = flag && (CGAL::FPU_get_rounding_mode() == CGAL::FPU_NEAREST);
   cout << "near   : " << (int) flag << endl;

   return (int) !flag;
}
