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
#include <CGAL/use.h>




template < typename AK>
void test_AK_1(){
  CGAL::IO::set_pretty_mode(std::cout);
  CGAL::IO::set_pretty_mode(std::cerr);

  typedef typename AK::Integer Integer;
  typedef typename AK::Rational Rational;
  CGAL_USE_TYPE(Rational);


  typedef CGAL::Polynomial<Integer> Poly;
  typedef CGAL::Polynomial_traits_d<Poly> PT;
  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Integer" << std::endl;
  std::cerr << "--------------------------------------" << std::endl;
  CGAL::Test_Pol::test_multiple_dimensions(PT());

}


template < typename AK>
void test_AK_2(){
  CGAL::IO::set_pretty_mode(std::cout);
  CGAL::IO::set_pretty_mode(std::cerr);

  typedef typename AK::Integer Integer;
  typedef typename AK::Rational Rational;
  CGAL_USE_TYPE(Integer);

  typedef CGAL::Polynomial<Rational> Poly;
  typedef CGAL::Polynomial_traits_d<Poly> PT;
  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Rational" << std::endl;
  std::cerr << "---------------------------------------" << std::endl;
  CGAL::Test_Pol::test_multiple_dimensions(PT());

}



template < typename AK>
void test_AK_4(){
  CGAL::IO::set_pretty_mode(std::cout);
  CGAL::IO::set_pretty_mode(std::cerr);

  typedef typename AK::Integer Integer;
  typedef typename AK::Rational Rational;


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


template < typename AK>
void test_AK_5(){
  CGAL::IO::set_pretty_mode(std::cout);
  CGAL::IO::set_pretty_mode(std::cerr);

  typedef typename AK::Integer Integer;
  typedef typename AK::Rational Rational;
  CGAL_USE_TYPE(Integer);


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


template < typename AK>
void test_AK_6(){
  CGAL::IO::set_pretty_mode(std::cout);
  CGAL::IO::set_pretty_mode(std::cerr);

  typedef typename AK::Integer Integer;
  typedef typename AK::Rational Rational;
  CGAL_USE_TYPE(Integer);
  CGAL_USE_TYPE(Rational);

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


int main(){
    // Set wrong rounding mode to test modular arithmetic
    CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);

#ifdef CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
    typedef CGAL::Arithmetic_kernel AK;
    test_AK_1<AK>();
    test_AK_2<AK>();
    test_AK_4<AK>();
    test_AK_5<AK>();
    test_AK_6<AK>();
#else
  std::cout << "No default arithmetic kernel has been found.\nNothing was tested" << std::endl;
#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

  return 0;
}
