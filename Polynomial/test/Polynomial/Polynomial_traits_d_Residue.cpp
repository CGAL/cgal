#include <iostream>
#include <cassert>

#include <CGAL/Arithmetic_kernel.h>

#ifdef CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

//#include <CGAL/FPU.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Test/_test_polynomial_traits_d.h>

int main()
{
  // Set wrong rounding mode to test modular arithmetic 
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);
  
  
  {
    //  Enforce IEEE double precision and to nearest before 
    //  using modular arithmetic
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);
    
    typedef CGAL::Polynomial< CGAL::Residue > Poly;
    typedef CGAL::Polynomial_traits_d<Poly> PT;
    std::cerr << std::endl;
    std::cerr << 
      "Test for coefficient type CGAL::Residue" 
              << std::endl;
    std::cerr << 
      "----------------------------------------------------------------------"
              << std::endl;    
    CGAL::Test_Pol::test_multiple_dimensions(PT());  
  }
  return 0;
} 
#else

int main()
{
  std::cout << "No default arithmetic kernel has been found.\nNothing was tested" << std::endl;
  return 0;
}
    
#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

