// This tests the rounding mode functions.

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>
#include <iostream>

typedef CGAL::Interval_nt_advanced NT_adv;
typedef CGAL::Interval_nt<>        NT;

void print_res (bool res)
{ std::cout << (res ? "ok" : "ERROR") << std::endl; }

// The results of 1-epsilon and -1+epsilon are enough
// to detect exactly the current rounding mode.
//                          1-MIN_DOUBLE
//                        +------+-------+
//                        |  1   | 1-ulp |
//               +--------+------+-------+
// -1+MIN_DOUBLE | -1     | near |  -inf |
//               | -1+ulp | +inf |  zero |
//               +--------+------+-------+

// I use a global variable here to avoid constant propagation.
double IA_min_double = CGAL_IA_MIN_DOUBLE;

CGAL::FPU_CW_t
FPU_empiric_test()
{
    double y = 1.0, z = -1.0;
    double ye, ze;
    ye = y - IA_min_double;
    ze = z + IA_min_double;
    if (y == ye && z == ze) return CGAL_FE_TONEAREST;
    if (y == ye) return CGAL_FE_UPWARD;
    if (z == ze) return CGAL_FE_DOWNWARD;
    return CGAL_FE_TOWARDZERO;
}

// This variable is global in order to stop constant propagation.
double m = 0.5;

CGAL::FPU_CW_t FPU_empiric_test_mul ()
{
    int i;
    for (i=0; i<10; i++) {m*=m; /* std::cout <<c << std::endl; */ }
    double a = m*m;
    double b = (-m)*m;
    std::cout << "m = " << m << "\n m*m = " << a << "\n (-m)*m = " << b;
    std::cout << std::endl;
// Note: it's not supposed to work here like that.
    if ((a == 0.0) && (b == 0.0)) return CGAL_FE_TONEAREST;
    if (a > 0.0) return CGAL_FE_UPWARD;
    if (b < 0.0) return CGAL_FE_DOWNWARD;
    return CGAL_FE_TOWARDZERO;
}

void print_rounding_name (CGAL::FPU_CW_t r)
{
  switch (r) {
  case CGAL_FE_TONEAREST:  std::cout << "NEAR\n"; break;
  case CGAL_FE_DOWNWARD:   std::cout << "DOWN\n"; break;
  case CGAL_FE_UPWARD:     std::cout << "UP\n"; break;
  case CGAL_FE_TOWARDZERO: std::cout << "ZERO\n"; break;
  default:                 std::cout << "unknown !\n";
  }
}

int main()
{
   bool flag = true;

   flag = flag && (FPU_empiric_test() == CGAL_FE_TONEAREST);
   std::cout << "default: ";
   print_res(flag);

   // Should be a no-op.
   CGAL::FPU_set_cw(CGAL::FPU_get_cw());
   flag = flag && (FPU_empiric_test() == CGAL_FE_TONEAREST);
   std::cout << "get/set: ";
   print_res(flag);
   if (!flag) print_rounding_name(FPU_empiric_test());

   // Rounding to zero.
   CGAL::FPU_set_cw(CGAL_FE_TOWARDZERO);
   flag = flag && (FPU_empiric_test() == CGAL_FE_TOWARDZERO);
   std::cout << "zero   : ";
   print_res(flag);
   if (!flag) print_rounding_name(FPU_empiric_test());

   // Rounding to infinity.
   CGAL::FPU_set_cw(CGAL_FE_UPWARD);
   flag = flag && (FPU_empiric_test() == CGAL_FE_UPWARD);
   std::cout << "+inf   : ";
   print_res(flag);
   if (!flag) print_rounding_name(FPU_empiric_test());

   // Rounding to minus infinity.
   CGAL::FPU_set_cw(CGAL_FE_DOWNWARD);
   flag = flag && (FPU_empiric_test() == CGAL_FE_DOWNWARD);
   std::cout << "-inf   : ";
   print_res(flag);
   if (!flag) print_rounding_name(FPU_empiric_test());

   // Rounding to nearest.
   CGAL::FPU_set_cw(CGAL_FE_TONEAREST);
   flag = flag && (FPU_empiric_test() == CGAL_FE_TONEAREST);
   std::cout << "near   : ";
   print_res(flag);
   if (!flag) print_rounding_name(FPU_empiric_test());

   return (int) !flag;
}
