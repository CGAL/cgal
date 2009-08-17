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





template < typename AT> 
void test_AT(){
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);

  typedef typename AT::Integer Integer;
  typedef typename AT::Rational Rational; 
  
  {
  typedef CGAL::Polynomial<Integer> Poly;
  typedef CGAL::Polynomial_traits_d<Poly> PT; 
  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Integer" << std::endl;
  std::cerr << "--------------------------------------" << std::endl;
  CGAL::Test_Pol::test_multiple_dimensions(PT());
  }
  {
  typedef CGAL::Polynomial<Rational> Poly;
  typedef CGAL::Polynomial_traits_d<Poly> PT; 
  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Rational" << std::endl;
  std::cerr << "---------------------------------------" << std::endl;
  CGAL::Test_Pol::test_multiple_dimensions(PT());
  }
  {
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
  {
  typedef CGAL::Polynomial< CGAL::Sqrt_extension< Rational, Integer > > Poly;
  typedef CGAL::Polynomial_traits_d<Poly> PT;
  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Sqrt_extension< Rational, Integer >"
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  CGAL::Test_Pol::test_multiple_dimensions(PT());    
  }
  {
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
  }
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
}    

int main(){

    // Set wrong rounding mode to test modular arithmetic 
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);

#ifdef CGAL_USE_LEDA
  {        
    typedef CGAL::LEDA_arithmetic_kernel AT;
    test_AT<AT>();
  }
#endif
#ifdef CGAL_USE_CORE
  {    
    typedef CGAL::CORE_arithmetic_kernel AT;
    test_AT<AT>();
  }
#endif

  return 0;
}
