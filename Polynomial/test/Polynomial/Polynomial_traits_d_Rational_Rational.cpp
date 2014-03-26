#include <iostream>
#include <cassert>

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>

#ifdef CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

#include <CGAL/Sqrt_extension.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Test/_test_polynomial_traits_d.h>

int main()
{
  typedef CGAL::Arithmetic_kernel AK; 
  typedef AK::Rational Rational; 
  
  typedef CGAL::Polynomial< CGAL::Sqrt_extension< Rational, Rational > > Poly;
  typedef CGAL::Polynomial_traits_d<Poly> PT;
  std::cerr << std::endl;
  std::cerr << 
    "Test for coefficient type Sqrt_extension< Rational, Rational >" 
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  CGAL::Test_Pol::test_multiple_dimensions(PT());    
  return 0;
}
 
#else
int main()
{
  std::cout << "No default arithmetic kernel has been found.\nNothing was tested" << std::endl;
  return 0;
}

#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
