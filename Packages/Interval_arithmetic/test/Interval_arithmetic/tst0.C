// This tests the rounding mode functions.

#include <CGAL/Interval_arithmetic.h>

using namespace CGAL;

// 5 temporary functions to test the inlining of the compiler.
bool triv_test_1 ()
{ return Interval_nt_advanced(1.0) < Interval_nt_advanced(2.0); }

bool triv_test_2 ()
{ return true; }

Interval_nt_advanced triv_mul_1 (Interval_nt_advanced x)
{ return x * Interval_nt_advanced(1.0); }

Interval_nt_advanced triv_1_mul (Interval_nt_advanced x)
{ return Interval_nt_advanced(1.0) * x; }

Interval_nt_advanced triv (Interval_nt_advanced x)
{ return x; }

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
    if (y == ye && z == ze) return FPU_cw_near;
    if (y == ye) return FPU_cw_up;
    if (z == ze) return FPU_cw_down;
    return FPU_cw_zero;
}

FPU_CW_t FPU_empiric_test_mul ()
{
    // If not marked "volatile", the result is false when optimizing
    // because the constants are pre-computed at compile time !!!
    // volatile const double m = CGAL_IA_MIN_DOUBLE;
    volatile double m = 0.5;
    int i;
    for (i=0; i<10; i++) {m*=m; /* cout <<c << endl; */ }
    double a = m*m;
    double b = (-m)*m;
    cout << "m = " << m << "\n m*m = " << a << "\n (-m)*m = " << b << endl;
// Note: it's not supposed to work here like that.
    if ((a == 0.0) && (b == 0.0)) return FPU_cw_near;
    if (a > 0.0) return FPU_cw_up;
    if (b < 0.0) return FPU_cw_down;
    return FPU_cw_zero;
}

void print_rounding_name (FPU_CW_t r)
{
  switch (r) {
  case FPU_cw_near: cout << "FPU_cw_near\n"; break;
  case FPU_cw_down: cout << "FPU_cw_down\n"; break;
  case FPU_cw_up:   cout << "FPU_cw_up\n"; break;
  case FPU_cw_zero: cout << "FPU_cw_zero\n"; break;
  default:          cout << "unknown !\n";
  }
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
   if (!flag) print_rounding_name(FPU_empiric_test());

   // Rounding to zero.
   FPU_set_cw(FPU_cw_zero);
   flag = flag && (FPU_empiric_test() == FPU_cw_zero);
   cout << "zero   : " << (int) flag << endl;
   if (!flag) print_rounding_name(FPU_empiric_test());

   // Rounding to infinity.
   FPU_set_cw(FPU_cw_up);
   flag = flag && (FPU_empiric_test() == FPU_cw_up);
   cout << "+inf   : " << (int) flag << endl;
   if (!flag) print_rounding_name(FPU_empiric_test());

   // Rounding to minus infinity.
   FPU_set_cw(FPU_cw_down);
   flag = flag && (FPU_empiric_test() == FPU_cw_down);
   cout << "-inf   : " << (int) flag << endl;
   if (!flag) print_rounding_name(FPU_empiric_test());

   // Rounding to nearest.
   FPU_set_cw(FPU_cw_near);
   flag = flag && (FPU_empiric_test() == FPU_cw_near);
   cout << "near   : " << (int) flag << endl;
   if (!flag) print_rounding_name(FPU_empiric_test());

   return (int) !flag;
}
