// This tests the rounding mode functions.

#include <CGAL/Interval_arithmetic.h>

typedef CGAL::Interval_nt_advanced NT_adv;

// 5 temporary functions to test the inlining of the compiler.
bool triv_test_1 ()
{ return NT_adv(1.0) < NT_adv(2.0); }

bool triv_test_2 ()
{ return true; }

NT_adv triv_mul_1 (NT_adv x)
{ return x * NT_adv(1.0); }

NT_adv triv_1_mul (NT_adv x)
{ return NT_adv(1.0) * x; }

NT_adv triv (NT_adv x)
{ return x; }

void print_res (bool res)
{ std::cout << (res ? "ok" : "ERROR") << std::endl; }

CGAL::FPU_CW_t FPU_empiric_test_mul ()
{
    // If not marked "volatile", the result is false when optimizing
    // because the constants are pre-computed at compile time !!!
    // volatile const double m = CGAL_IA_MIN_DOUBLE;
    volatile double m = 0.5;
    int i;
    for (i=0; i<10; i++) {m*=m; /* std::cout <<c << std::endl; */ }
    double a = m*m;
    double b = (-m)*m;
    std::cout << "m = " << m << "\n m*m = " << a << "\n (-m)*m = " << b;
    std::cout << std::endl;
// Note: it's not supposed to work here like that.
    if ((a == 0.0) && (b == 0.0)) return CGAL::FPU_cw_near;
    if (a > 0.0) return CGAL::FPU_cw_up;
    if (b < 0.0) return CGAL::FPU_cw_down;
    return CGAL::FPU_cw_zero;
}

void print_rounding_name (CGAL::FPU_CW_t r)
{
  switch (r) {
  case CGAL::FPU_cw_near: std::cout << "FPU_cw_near\n"; break;
  case CGAL::FPU_cw_down: std::cout << "FPU_cw_down\n"; break;
  case CGAL::FPU_cw_up:   std::cout << "FPU_cw_up\n"; break;
  case CGAL::FPU_cw_zero: std::cout << "FPU_cw_zero\n"; break;
  default:          std::cout << "unknown !\n";
  }
}

int main()
{
   bool flag = true;

   flag = flag && (CGAL::FPU_empiric_test() == CGAL::FPU_cw_near);
   std::cout << "default: ";
   print_res(flag);

   // Should be a no-op.
   CGAL::FPU_set_cw(CGAL::FPU_get_cw());
   flag = flag && (CGAL::FPU_empiric_test() == CGAL::FPU_cw_near);
   std::cout << "get/set: ";
   print_res(flag);
   if (!flag) print_rounding_name(CGAL::FPU_empiric_test());

   // Rounding to zero.
   CGAL::FPU_set_cw(CGAL::FPU_cw_zero);
   flag = flag && (CGAL::FPU_empiric_test() == CGAL::FPU_cw_zero);
   std::cout << "zero   : ";
   print_res(flag);
   if (!flag) print_rounding_name(CGAL::FPU_empiric_test());

   // Rounding to infinity.
   CGAL::FPU_set_cw(CGAL::FPU_cw_up);
   flag = flag && (CGAL::FPU_empiric_test() == CGAL::FPU_cw_up);
   std::cout << "+inf   : ";
   print_res(flag);
   if (!flag) print_rounding_name(CGAL::FPU_empiric_test());

   // Rounding to minus infinity.
   CGAL::FPU_set_cw(CGAL::FPU_cw_down);
   flag = flag && (CGAL::FPU_empiric_test() == CGAL::FPU_cw_down);
   std::cout << "-inf   : ";
   print_res(flag);
   if (!flag) print_rounding_name(CGAL::FPU_empiric_test());

   // Rounding to nearest.
   CGAL::FPU_set_cw(CGAL::FPU_cw_near);
   flag = flag && (CGAL::FPU_empiric_test() == CGAL::FPU_cw_near);
   std::cout << "near   : ";
   print_res(flag);
   if (!flag) print_rounding_name(CGAL::FPU_empiric_test());

   return (int) !flag;
}
