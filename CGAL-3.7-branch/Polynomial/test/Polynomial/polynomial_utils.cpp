// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : test/Polynomial/polynomial_utils.cpp
// CGAL_release   : $Name:  $
// Revision      : $Revision: 46395 $
// Revision_date : $Date: 2008-10-21 14:59:59 +0200 (Tue, 21 Oct 2008) $
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

#include <CGAL/polynomial_utils.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/Exponent_vector.h>
#include <CGAL/Sqrt_extension.h>


template <typename AK>
void test_polynomial_utils(){
  //CGAL::set_pretty_mode(std::cout);
  
  typedef typename AK::Integer Integer;
  typedef typename AK::Rational Rational; 
//  typedef CGAL::Sqrt_extension<Integer,Integer> EXT_INT;
//  typedef CGAL::Sqrt_extension<Rational,Integer> EXT_RAT;

  typedef typename CGAL::Polynomial_type_generator<Integer,3>::Type POLY_INT_3;
  typedef CGAL::Polynomial_traits_d<POLY_INT_3>     PT_3;
  typedef typename PT_3::Innermost_coefficient_type ICOEFF;
  typedef typename PT_3::Coefficient_type           COEFF;
  
  typename CGAL::Polynomial_traits_d<POLY_INT_3>::Shift shift; 
  POLY_INT_3 x = shift(POLY_INT_3(1),1,0);
  POLY_INT_3 y = shift(POLY_INT_3(1),1,1);
  POLY_INT_3 z = shift(POLY_INT_3(1),1,2);

  POLY_INT_3 p = -5*x*x*x*y+7*z*z*y;  
  //std::cout << p << std::endl; 

// GetCoefficient
  assert(CGAL::get_coefficient(x*x*y+z,1)== 1);
  assert(CGAL::get_coefficient(x*x*y+z,0)== x*x*y);
// GetInnermostCoefficient
  assert(CGAL::get_innermost_coefficient(x*x*y+z,CGAL::Exponent_vector(0,0,1))== 1);
  assert(CGAL::get_innermost_coefficient(x*x*y+z,CGAL::Exponent_vector(2,1,0))== 1);
  assert(CGAL::get_innermost_coefficient(x*x*y+z,CGAL::Exponent_vector(2,1,1))== 0);
// ConstructCoefficientConstIteratorRange
// ConstructInnermostCoefficientConstIteratorRange
// Swap
  assert(CGAL::swap(x*x*y,0,1)==x*y*y);
// Move
  assert(CGAL::move(x*x*y,0,2)==x*z*z);
// Permute
  int permutation[3] = {1,2,0};
  assert(CGAL::permute(x*x*y,permutation,permutation+3)==y*y*z);
// Degree
  assert(CGAL::degree(p) == 2);
  assert(CGAL::degree(p,0) == 3);
// TotalDegree
  assert(CGAL::total_degree(p)  == 4);
// DegreeVector
  assert(CGAL::degree_vector(p) == CGAL::Exponent_vector(0,1,2));
// LeadingCoefficient
  assert(CGAL::leading_coefficient(p)   == 7*y);
// InnermostLeadingCoefficient
 assert(CGAL::innermost_leading_coefficient(p) == 7);
// Canonicalize
  assert(CGAL::canonicalize(2*p) == p);
  assert(CGAL::is_zero(CGAL::canonicalize(POLY_INT_3(0))));
// Differentiate
  assert(CGAL::differentiate(p)   ==  14*y*z);
  assert(CGAL::differentiate(p,0) ==  -15*x*x*y);
  assert(CGAL::differentiate(p,1) ==  7*z*z-5*x*x*x);
  assert(CGAL::differentiate(p,2) ==  14*y*z);
// Evaluate
  assert(CGAL::evaluate(p,COEFF(2)) ==  (-5*x*x*x + 28)*y );
// EvaluateHomogeneous
  assert(CGAL::evaluate_homogeneous(p,COEFF(2),COEFF(3)) 
      == (-45*x*x*x + 28)*y );
// Substitute
  std::vector<Rational> vec; 
  vec.push_back(Rational(1));
  vec.push_back(Rational(2));
  vec.push_back(Rational(3));
  assert(CGAL::substitute(p,vec.begin(), vec.end()) == Rational(116)); 
// IsZeroAt
  assert(CGAL::is_zero_at(p,vec.begin(), vec.end()) == false); 
// SignAt
  assert(CGAL::sign_at(p,vec.begin(), vec.end()) == CGAL::POSITIVE);
  vec.push_back(Rational(4));
// SubstituteHomogeneous
  assert(CGAL::substitute_homogeneous(p,vec.begin(), vec.end()) == Rational(494)); 
// IsZeroAtHomogeneous
  assert(CGAL::is_zero_at_homogeneous(p,vec.begin(), vec.end()) == false); 
// SignAtHomogeneous
  assert(CGAL::sign_at_homogeneous(p,vec.begin(), vec.end()) == CGAL::POSITIVE);

// Compare
  assert(CGAL::compare(p,p)  == CGAL::EQUAL);
  assert(CGAL::compare(p,-p) == CGAL::LARGER);
  assert(CGAL::compare(p,2*p) == CGAL::SMALLER);
  
// UnivariateContent
  assert(CGAL::univariate_content(p) == y);
// MultivariateContent
  assert(CGAL::multivariate_content(2*p) == 2);
  assert(CGAL::multivariate_content(-12*p) == 12);

// SquareFreeFactorize
  {
    std::vector<std::pair<POLY_INT_3,int> > sqff_vec;
    CGAL::square_free_factorize(p*p,std::back_inserter(sqff_vec));
    POLY_INT_3 tmp(1);
    assert(sqff_vec.size() >= 2);  
    for(unsigned int i = 0; i < sqff_vec.size();i++){
      tmp *= CGAL::ipower(sqff_vec[i].first,sqff_vec[i].second);
    }
    assert(tmp == p*p);
  }
// MakeSquareFree
  assert(CGAL::make_square_free(p*p*y) == p);
// IsSquareFree
  assert(CGAL::is_square_free(p*p*y) == false);
  assert(CGAL::is_square_free(p) == true);
  
  
// PseudoDivision
  {
    POLY_INT_3 q,r;
    COEFF D;
    POLY_INT_3 g = 5*z-x*y*z;
    CGAL::pseudo_division(p,g,q,r,D);
    assert(D*p == q*g+r);  
// PseudoDivisionQuotient
    assert(CGAL::pseudo_division_quotient(p,g) == q); 
// PseudoDivisionRemainder
    assert(CGAL::pseudo_division_remainder(p,g) == r); 
  }
// GcdUpToConstantFactor
  assert(CGAL::gcd_up_to_constant_factor(5*p,5*y) == y);
// IntegralDivisionUpToConstantFactor
  assert(CGAL::integral_division_up_to_constant_factor(-5*p,y) 
      == CGAL::integral_division(p,y));
// UnivariateContentUpToConstantFactor
  assert(CGAL::univariate_content_up_to_constant_factor(-5*p) == y);
// SquareFreeFactorizeUpToConstantFactor
  {
    std::vector<std::pair<POLY_INT_3,int> > sqff_vec;
    CGAL::square_free_factorize_up_to_constant_factor
      (25*p*p,std::back_inserter(sqff_vec));
    POLY_INT_3 tmp(1);
    assert(sqff_vec.size() >= 2);  
    for(unsigned int i = 0; i < sqff_vec.size();i++){
      tmp *= CGAL::ipower(sqff_vec[i].first,sqff_vec[i].second);
    }
    assert(CGAL::canonicalize(tmp) == CGAL::canonicalize(p*p));
  }

//Shift
  assert(x == CGAL::shift(POLY_INT_3(1),1,0));
  assert(y == CGAL::shift(POLY_INT_3(1),1,1));
  assert(z == CGAL::shift(POLY_INT_3(1),1,2));
  assert(z == CGAL::shift(POLY_INT_3(1),1));

//Negate
  // p = -5*x*x*x*y+7*z*z*y
  assert(CGAL::negate(p,0) ==  5*x*x*x*y+7*z*z*y);
  assert(CGAL::negate(p,1) ==  5*x*x*x*y-7*z*z*y);
  assert(CGAL::negate(p,2) == -5*x*x*x*y+7*z*z*y);
  assert(CGAL::negate(p)   == -5*x*x*x*y+7*z*z*y);
//Invert
  assert(CGAL::invert(p,0) == -5*y+7*z*z*y*x*x*x);
  assert(CGAL::invert(p,1) == -5*x*x*x+7*z*z);
  assert(CGAL::invert(p,2) == -5*x*x*x*y*z*z+7*y);
  assert(CGAL::invert(p)   == -5*x*x*x*y*z*z+7*y);
//Translate
  assert(CGAL::translate(x*y*z,ICOEFF(2),0) == (x+2)*y*z);
  assert(CGAL::translate(x*y*z,ICOEFF(2),1) == (y+2)*x*z);
  assert(CGAL::translate(x*y*z,ICOEFF(2),2) == (z+2)*x*y);
  assert(CGAL::translate(x*y*z,ICOEFF(2)) == (z+2)*x*y);
//TranslateHomogeneous
  assert(CGAL::translate_homogeneous(x*y*z,ICOEFF(2),ICOEFF(3),0)==(3*x+2)*y*z);
  assert(CGAL::translate_homogeneous(x*y*z,ICOEFF(2),ICOEFF(3),1)==(3*y+2)*x*z);
  assert(CGAL::translate_homogeneous(x*y*z,ICOEFF(2),ICOEFF(3),2)==(3*z+2)*x*y);
  assert(CGAL::translate_homogeneous(x*y*z,ICOEFF(2),ICOEFF(3))  ==(3*z+2)*x*y);
//Scale
  assert(CGAL::scale(x*x+y*y+z*z,ICOEFF(2),0) == 4*x*x+y*y+z*z);
  assert(CGAL::scale(x*x+y*y+z*z,ICOEFF(2),1) == x*x+4*y*y+z*z);
  assert(CGAL::scale(x*x+y*y+z*z,ICOEFF(2),2) == x*x+y*y+4*z*z);
  assert(CGAL::scale(x*x+y*y+z*z,ICOEFF(2)) ==   x*x+y*y+4*z*z);
//ScaleHomogeneous
  assert(CGAL::scale_homogeneous(x*x+y*y+z*z,ICOEFF(2),ICOEFF(3),0) 
      == 4*x*x+9*y*y+9*z*z);
  assert(CGAL::scale_homogeneous(x*x+y*y+z*z,ICOEFF(2),ICOEFF(3),1) 
      == 9*x*x+4*y*y+9*z*z);
  assert(CGAL::scale_homogeneous(x*x+y*y+z*z,ICOEFF(2),ICOEFF(3),2) 
      == 9*x*x+9*y*y+4*z*z);
  assert(CGAL::scale_homogeneous(x*x+y*y+z*z,ICOEFF(2),ICOEFF(3)) 
      == 9*x*x+9*y*y+4*z*z);
//Resultant
  assert(CGAL::is_zero(CGAL::resultant(p,p)));
  assert(CGAL::resultant(3*x*x*z+x*y,5*y*y*z+y*x)
      == -y*x*(5*y*y-3*x*x)); // Maple ;-)
}

template <typename AK>
void test_canonicalize(){
  
  typedef typename AK::Integer Integer;
  typedef typename AK::Rational Rational; 
  typedef CGAL::Sqrt_extension<Integer,Integer> EXT_INT;
  typedef CGAL::Sqrt_extension<Rational,Integer> EXT_RAT;

  {
    typedef Integer NT; 
    typedef typename CGAL::Polynomial_type_generator<NT,2>::Type POLY_2;
    POLY_2 x = CGAL::shift(POLY_2(1),1,0);
    POLY_2 y = CGAL::shift(POLY_2(1),1,1);
    POLY_2 p = -5*x*x*x*y+7*y;
    POLY_2 q = CGAL::canonicalize(2 * p); 
    assert(CGAL::innermost_leading_coefficient(q) == NT(5));
    assert(CGAL::canonicalize(p) == CGAL::canonicalize(-p));
  }{
    typedef Rational NT; 
    typedef typename CGAL::Polynomial_type_generator<NT,2>::Type POLY_2;
    POLY_2 x = CGAL::shift(POLY_2(1),1,0);
    POLY_2 y = CGAL::shift(POLY_2(1),1,1);
    POLY_2 p = -5*x*x*x*y+7*y;
    POLY_2 q = CGAL::canonicalize(2 * p); 
    assert(CGAL::innermost_leading_coefficient(q) == NT(1));
    assert(CGAL::canonicalize(p) == CGAL::canonicalize(-p));
  }{
    typedef EXT_INT NT; 
    typedef typename CGAL::Polynomial_type_generator<NT,2>::Type POLY_2;
    POLY_2   x = CGAL::shift(POLY_2(1),1,0);
    POLY_2   y = CGAL::shift(POLY_2(1),1,1);
    POLY_2   p = -5*x*x*x*y+7*y;
    EXT_INT ex(Integer(2),Integer(5),Integer(7));
    POLY_2   q = CGAL::canonicalize(ex * p); 
    assert(CGAL::innermost_leading_coefficient(q) == NT(5));
    assert(CGAL::canonicalize(p) == CGAL::canonicalize(-p));
  }{
    typedef EXT_RAT NT; 
    typedef typename CGAL::Polynomial_type_generator<NT,2>::Type POLY_2;
    POLY_2   x = CGAL::shift(POLY_2(1),1,0);
    POLY_2   y = CGAL::shift(POLY_2(1),1,1);
    POLY_2   p = -5*x*x*x*y+7*y;
    EXT_RAT ex(Integer(2),Integer(5),Integer(7));
    POLY_2   q = CGAL::canonicalize(ex * p); 
    assert(CGAL::innermost_leading_coefficient(q) == NT(1));
    assert(CGAL::canonicalize(p) == CGAL::canonicalize(-p));
  }
}


int main(){
#if CGAL_HAS_DEFAULT_ARITHMETIC_KERNEL 
  typedef CGAL::Arithmetic_kernel AK; 
  test_polynomial_utils<AK>();
  test_canonicalize<AK>();
#endif 

}
