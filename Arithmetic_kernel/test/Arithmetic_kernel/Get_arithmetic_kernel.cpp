
#include <iostream>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Get_arithmetic_kernel.h>
#include <CGAL/use.h>

#if defined(CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL)

#include <CGAL/Test/_test_arithmetic_kernel.h>

int main() {

  typedef CGAL::Arithmetic_kernel AK;
  typedef AK::Integer Integer;
  typedef AK::Rational Rational;
  typedef AK::Bigfloat_interval BFI;
  {
    typedef CGAL::Get_arithmetic_kernel<Integer>::Arithmetic_kernel AK_;
    CGAL_USE_TYPE(AK_);
    static_assert(std::is_same<AK,AK_>::value);
  }
  {
    typedef CGAL::Get_arithmetic_kernel<Rational>::Arithmetic_kernel AK_;
    CGAL_USE_TYPE(AK_);
    static_assert(std::is_same<AK,AK_>::value);
  }
  {
    typedef CGAL::Get_arithmetic_kernel<BFI>::Arithmetic_kernel AK_;
    CGAL_USE_TYPE(AK_);
    static_assert(std::is_same<AK,AK_>::value);
  }
  return 0;
}

#else
#warning CGAL has no default CGAL::Arithmetic kernel
int main() { return 0; }
#endif
