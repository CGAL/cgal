// This tests the rounding mode functions.

#include <CGAL/Interval_arithmetic.h>

using namespace CGAL;

// Rounding mode empiric testing.

// The results of 1-epsilon and -1+epsilon are enough
// to detect exactly the rounding mode.
// (epsilon = MIN_DOUBLE, ulp = 2^-52 or 2^-53).
// ----------------------------------------------------
// rounding mode:        +inf    -inf    0       nearest
// ----------------------------------------------------
//  1-epsilon            1       1-ulp   1-ulp   1
// -1+epsilon           -1+ulp  -1      -1+ulp  -1
// ----------------------------------------------------

FPU_CW_t FPU_empiric_test ()
{
    // If not marked "volatile", the result is false when optimizing
    // because the constants are pre-computed at compile time !!!
    volatile const double m = CGAL_IA_MIN_DOUBLE;
    const double y = 1.0, z = -1.0;
    double ye, ze;
    ye = y - m;
    ze = z + m;
    if ((y == ye) && (z == ze)) return FPU_cw_near;
    if (y == ye) return FPU_cw_up;
    if (z == ze) return FPU_cw_down;
    return FPU_cw_zero;
}

int main()
{
   bool flag = true;

   flag = flag && (FPU_empiric_test() == FPU_cw_near);
   cout << "default: " << (int) flag << endl;

   // Should be a no-op.
   FPU_set_cw(FPU_get_cw());
   flag = flag && (FPU_empiric_test() == FPU_cw_near);
   cout << "get/set: " << (int) flag << endl;

   // Rounding to zero.
   FPU_set_cw(FPU_cw_zero);
   flag = flag && (FPU_empiric_test() == FPU_cw_zero);
   cout << "zero   : " << (int) flag << endl;

   // Rounding to infinity.
   FPU_set_cw(FPU_cw_up);
   flag = flag && (FPU_empiric_test() == FPU_cw_up);
   cout << "+inf   : " << (int) flag << endl;

   // Rounding to minus infinity.
   FPU_set_cw(FPU_cw_down);
   flag = flag && (FPU_empiric_test() == FPU_cw_down);
   cout << "-inf   : " << (int) flag << endl;

   // Rounding to nearest.
   FPU_set_cw(FPU_cw_near);
   flag = flag && (FPU_empiric_test() == FPU_cw_near);
   cout << "near   : " << (int) flag << endl;

   return (int) !flag;
}
