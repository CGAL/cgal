#include <CGAL/config.h>
#include <iostream>

#ifdef CGAL_USE_GMP

#include <CGAL/Gmpz.h>
#include <CGAL/Polynomial.h>

// Function in case  Polynomial is Modularizable
template< typename Polynomial >
bool may_have_common_factor(
    const Polynomial& p1, const Polynomial& p2, CGAL::Tag_true){
  std::cout<< "The type is modularizable" << std::endl;

  // Enforce IEEE double precision and rounding mode to nearest
  // before useing modular arithmetic
  CGAL::Protect_FPU_rounding<true> pfr(CGAL_FE_TONEAREST);

  // Use Modular_traits to convert to polynomials with modular coefficients
  typedef CGAL::Modular_traits<Polynomial> MT;
  typedef typename MT::Residue_type MPolynomial;
  typedef typename MT::Modular_image Modular_image;
  MPolynomial mp1 = Modular_image()(p1);
  MPolynomial mp2 = Modular_image()(p2);

  // check for unlucky primes, the polynomials should not lose a degree
  typename CGAL::Polynomial_traits_d<Polynomial>::Degree  degree;
  typename CGAL::Polynomial_traits_d<MPolynomial>::Degree mdegree;
  if ( degree(p1) != mdegree(mp1)) return true;
  if ( degree(p2) != mdegree(mp2)) return true;

  // compute gcd for modular images
  MPolynomial mg  = CGAL::gcd(mp1,mp2);

  // if the modular gcd is not trivial: return true
  if ( mdegree(mg) > 0 ){
    std::cout << "The gcd may be non trivial" << std::endl;
    return true;
  }else{
    std::cout << "The gcd is trivial" << std::endl;
    return false;
  }
}

// This function returns true, since the filter is not applicable
template< typename Polynomial >
bool may_have_common_factor(
    const Polynomial&, const Polynomial&, CGAL::Tag_false){
  std::cout<< "The type is not modularizable" << std::endl;
  return true;
}

template< typename Polynomial >
Polynomial modular_filtered_gcd(const Polynomial& p1, const Polynomial& p2){
  typedef CGAL::Modular_traits<Polynomial> MT;
  typedef typename MT::Is_modularizable Is_modularizable;

  // Try to avoid actual gcd computation
  if (may_have_common_factor(p1,p2, Is_modularizable())){
    // Compute gcd, since the filter indicates a common factor
    return CGAL::gcd(p1,p2);
  }else{
    typename CGAL::Polynomial_traits_d<Polynomial>::Univariate_content  content;
    typename CGAL::Polynomial_traits_d<Polynomial>::Construct_polynomial construct;
    return construct(CGAL::gcd(content(p1),content(p2))); // return trivial gcd
  }
}

int main(){
  CGAL::set_pretty_mode(std::cout);

  typedef CGAL::Gmpz NT;
  typedef CGAL::Polynomial<NT> Poly;
  CGAL::Polynomial_traits_d<Poly>::Construct_polynomial construct;

  Poly  f1=construct(NT(2), NT(6), NT(4));
  Poly  f2=construct(NT(12), NT(4), NT(8));
  Poly  f3=construct(NT(3), NT(4));

  std::cout << "f1        : " << f1 << std::endl;
  std::cout << "f2        : " << f2 << std::endl;

  std::cout << "compute modular filtered gcd(f1,f2): " << std::endl;
  Poly g1 = modular_filtered_gcd(f1,f2);
  std::cout << "gcd(f1,f2): " << g1 << std::endl;

  std::cout << std::endl;
  Poly p1 = f1*f3;
  Poly p2 = f2*f3;

  std::cout << "f3        : " << f3 << std::endl;
  std::cout << "p1=f1*f3  : " << p1 << std::endl;
  std::cout << "p2=f2*f3  : " << p2 << std::endl;

  std::cout << "compute modular filtered gcd(p1,p2): " << std::endl;
  Poly g2 = modular_filtered_gcd(p1,p2);
  std::cout << "gcd(p1,p2): " << g2 << std::endl;
}

#else

int main (){
  std::cout << " This example needs GMP! " << std::endl;
}

#endif
