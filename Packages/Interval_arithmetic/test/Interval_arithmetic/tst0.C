// This tests the rounding mode.

#include <CGAL/Interval_arithmetic.h>

using namespace CGAL;

int main()
{
   bool flag = true;

   flag = flag && (FPU_get_rounding_mode() == FPU_NEAREST);
   cout << "default: " << (int) flag << endl;

   // FPU_set_rounding_to_zero();
   FPU_set_rounding_mode(FPU_ZERO);
   flag = flag && (FPU_get_rounding_mode() == FPU_ZERO);
   cout << "zero   : " << (int) flag << endl;

   // FPU_set_rounding_to_infinity();
   FPU_set_rounding_mode(FPU_PLUS_INFINITY);
   flag = flag && (FPU_get_rounding_mode() == FPU_PLUS_INFINITY);
   cout << "+inf   : " << (int) flag << endl;

   // FPU_set_rounding_to_minus_infinity();
   FPU_set_rounding_mode(FPU_MINUS_INFINITY);
   flag = flag && (FPU_get_rounding_mode() == FPU_MINUS_INFINITY);
   cout << "-inf   : " << (int) flag << endl;

   // FPU_set_rounding_to_nearest();
   FPU_set_rounding_mode(FPU_NEAREST);
   flag = flag && (FPU_get_rounding_mode() == FPU_NEAREST);
   cout << "near   : " << (int) flag << endl;

   return (int) !flag;
}
