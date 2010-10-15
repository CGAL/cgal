
// currently using modular arithmetic is far to slow,
// that is, we use interpolation for multivariate resultants 
// but no modular arithmetic. 

// These are the defaults 
//#define CGAL_RESULTANT_USE_MODULAR_ARITHMETIC 0
//#define CGAL_RESULTANT_USE_DECOMPOSE 1

#include <CGAL/basic.h>

#include <cassert>
#include <iostream>

#include <CGAL/Polynomial/resultant.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Residue.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Random.h>
#include <CGAL/Timer.h>

#include <CGAL/gen_sparse_polynomial.h>

static CGAL::Random my_rnd(346); // some seed 

template<class Polynomial_d> 
void test_resultant(){
  typedef typename CGAL::Polynomial_traits_d<Polynomial_d> PT;
  typedef typename PT::Coefficient_type Coeff; 
  typename PT::Resultant resultant;
   
  {
    Polynomial_d A(0);
    Polynomial_d B(0);
    assert(resultant(A,B)==Coeff(0));
  }{
    Polynomial_d A(4);
    Polynomial_d B(8);
    assert(resultant(A,B)==Coeff(1));
  }{
    Polynomial_d A(4);
    Polynomial_d B(Coeff(3), Coeff(5), Coeff(7));
    assert(resultant(A,B)==Coeff(16)); // lcoeff(A)^degree(B)
  }{
    Polynomial_d f(Coeff(2),Coeff(7),Coeff(1),Coeff(8),Coeff(1),Coeff(8));
    Polynomial_d g(Coeff(3),Coeff(1),Coeff(4),Coeff(1),Coeff(5),Coeff(9));
    assert(resultant(f,g) == Coeff(230664271L)); // Maple
        
    Polynomial_d h(Coeff(3),Coeff(4),Coeff(7),Coeff(7));
    Polynomial_d fh(f*h);
    Polynomial_d gh(g*h);
    assert(resultant(fh,gh) == Coeff(0) );
  } 
  for (int k = 0; k < 1; k++){
    Polynomial_d F2 = 
      CGAL::generate_sparse_random_polynomial<Polynomial_d>(my_rnd,5); 
    Polynomial_d F1 = 
      CGAL::generate_sparse_random_polynomial<Polynomial_d>(my_rnd,5);
    Coeff F = resultant(F1,F2); (void) F;
  }
  //std::cout << "end resultant: " << PT::d << std::endl;
}
   
#if CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL

int main(){

  // Set wrong rounding mode to test modular arithmetic 
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_UPWARD);

  CGAL::set_pretty_mode(std::cout);
  
  typedef CGAL::Arithmetic_kernel AK;
  typedef AK::Integer Integer; 
  typedef CGAL::Polynomial<Integer>      Polynomial_1;
  typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;
  typedef CGAL::Polynomial<Polynomial_2> Polynomial_3;

  test_resultant<Polynomial_1>();
  test_resultant<Polynomial_2>();
  test_resultant<Polynomial_3>();
 
  typedef CGAL::Polynomial<AK::Rational>   RPolynomial_1;
  typedef CGAL::Polynomial<RPolynomial_1>  RPolynomial_2;
  typedef CGAL::Polynomial<RPolynomial_2>  RPolynomial_3;

  test_resultant<RPolynomial_1>();
  test_resultant<RPolynomial_2>();
  test_resultant<RPolynomial_3>();  
    
  typedef CGAL::Sqrt_extension<AK::Integer, AK::Integer> EXT;
  typedef CGAL::Polynomial<EXT>            EPolynomial_1;
  typedef CGAL::Polynomial<EPolynomial_1>  EPolynomial_2;
  typedef CGAL::Polynomial<EPolynomial_2>  EPolynomial_3;

  test_resultant<EPolynomial_1>();
  test_resultant<EPolynomial_2>();
  test_resultant<EPolynomial_3>();
  
  typedef CGAL::Sqrt_extension<EXT, AK::Integer> EXT2;
  typedef CGAL::Polynomial<EXT2>            E2Polynomial_1;
  test_resultant<E2Polynomial_1>();


  // Enforce IEEE double precision and to nearest before using modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr2(CGAL_FE_TONEAREST);

  typedef CGAL::Polynomial<CGAL::Residue>  MPolynomial_1;
  typedef CGAL::Polynomial<MPolynomial_1>  MPolynomial_2;
  typedef CGAL::Polynomial<MPolynomial_2>  MPolynomial_3;

  test_resultant<MPolynomial_1>();
  test_resultant<MPolynomial_2>();
  test_resultant<MPolynomial_3>();
}

#else

int main(){
  std::cout << " Test needs a default arithmetic kernel " << std::endl;  
  return 0; 
}

#endif // CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL
 

