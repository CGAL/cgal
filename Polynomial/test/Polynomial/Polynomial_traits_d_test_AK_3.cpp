#include <iostream>
#include <cassert>

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/ipower.h>
#include <CGAL/Random.h>
#include <cmath>
#include <CGAL/Test/_test_polynomial_traits_d.h>





template < typename AK> 
void test_AK_3(){
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);

  typedef typename AK::Integer Integer;
  typedef typename AK::Rational Rational; 
  

  typedef CGAL::Polynomial< CGAL::Sqrt_extension< Integer, Integer > > Poly;
  typedef CGAL::Polynomial_traits_d<Poly> PT;  
  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Sqrt_extension< Integer, Integer >" 
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  CGAL::Test_Pol::test_multiple_dimensions(PT());    


}


int main(){
    // Set wrong rounding mode to test modular arithmetic 
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);

#ifdef CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
    typedef CGAL::Arithmetic_kernel AK; 
    test_AK_3<AK>();
#else
  std::cout << "No default arithmetic kernel has been found.\nNothing was tested" << std::endl;
#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

  return 0;
}
