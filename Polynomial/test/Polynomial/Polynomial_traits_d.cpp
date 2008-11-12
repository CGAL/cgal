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
#include <CGAL/Test/_test_Polynomial_traits_d.h>




template < typename AK>
void test_rebind(){
  typedef typename AK::Integer Integer; 
  typedef typename AK::Rational Rational;
  typedef CGAL::Polynomial<Integer> Poly_int_1;                
  typedef CGAL::Polynomial<Poly_int_1> Poly_int_2;              
  typedef CGAL::Polynomial<Poly_int_2> Poly_int_3;              
  typedef CGAL::Polynomial<Rational> Poly_rat_1;               
  typedef CGAL::Polynomial<Poly_rat_1> Poly_rat_2;              
  typedef CGAL::Polynomial<Poly_rat_2> Poly_rat_3;         
    
  typedef CGAL::Polynomial_traits_d<Poly_int_1> PT_int_1;
  typedef CGAL::Polynomial_traits_d<Poly_int_2> PT_int_2;
  typedef CGAL::Polynomial_traits_d<Poly_int_3> PT_int_3;
  typedef CGAL::Polynomial_traits_d<Poly_rat_1> PT_rat_1;
  typedef CGAL::Polynomial_traits_d<Poly_rat_2> PT_rat_2;
  typedef CGAL::Polynomial_traits_d<Poly_rat_3> PT_rat_3;

  typedef typename PT_int_1:: template Rebind<Integer,1>::Other PT_int_1_;
  typedef typename PT_int_3:: template Rebind<Integer,2>::Other PT_int_2_;
  typedef typename PT_rat_3:: template Rebind<Integer,3>::Other PT_int_3_;
  typedef typename PT_int_1:: template Rebind<Rational,1>::Other PT_rat_1_;
  typedef typename PT_rat_2:: template Rebind<Rational,2>::Other PT_rat_2_;
  typedef typename PT_int_2:: template Rebind<Rational,3>::Other PT_rat_3_;
    
  BOOST_STATIC_ASSERT((boost::is_same<PT_int_1_,PT_int_1>::value));
  BOOST_STATIC_ASSERT((boost::is_same<PT_int_2_,PT_int_2>::value));
  BOOST_STATIC_ASSERT((boost::is_same<PT_int_3_,PT_int_3>::value));
  BOOST_STATIC_ASSERT((boost::is_same<PT_rat_1_,PT_rat_1>::value));
  BOOST_STATIC_ASSERT((boost::is_same<PT_rat_2_,PT_rat_2>::value));
  BOOST_STATIC_ASSERT((boost::is_same<PT_rat_3_,PT_rat_3>::value));

  BOOST_STATIC_ASSERT((!boost::is_same<PT_rat_3_,PT_rat_2>::value));
}


template < typename AT> 
void test_AT(){
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);

  typedef typename AT::Integer Integer;
  typedef typename AT::Rational Rational; 

  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Integer" << std::endl;
  std::cerr << "--------------------------------------" << std::endl;
  CGAL::Test_Pol::test_multiple_dimensions<Integer>();

  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Rational" << std::endl;
  std::cerr << "---------------------------------------" << std::endl;
  CGAL::Test_Pol::test_multiple_dimensions<Rational>();
    
  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Sqrt_extension< Integer, Integer >" 
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  CGAL::Test_Pol::test_multiple_dimensions< CGAL::Sqrt_extension< Integer, Integer > >();    

  std::cerr << std::endl;
  std::cerr << "Test for coefficient type Sqrt_extension< Rational, Integer >"
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  CGAL::Test_Pol::test_multiple_dimensions< CGAL::Sqrt_extension< Rational, Integer > >();    

  std::cerr << std::endl;
  std::cerr << 
    "Test for coefficient type Sqrt_extension< Rational, Rational >" 
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
  CGAL::Test_Pol::test_multiple_dimensions< CGAL::Sqrt_extension< Rational, Rational > >();   

  test_rebind<AT>();

  std::cerr << std::endl;
  std::cerr << 
    "Test for coefficient type CGAL::Residue" 
            << std::endl;
  std::cerr << 
    "----------------------------------------------------------------------"
            << std::endl;    
 //  Enforce IEEE double precision before using modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);
  CGAL::Test_Pol::test_multiple_dimensions< CGAL::Residue >();   
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
