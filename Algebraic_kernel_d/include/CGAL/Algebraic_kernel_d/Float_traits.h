// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : 
//
// ============================================================================

// TODO: Some comments are original EXACUS comments and aren't adapted. So
//         they may be wrong now.

// TODO: should exponent type be long or Integer ? 

#ifndef CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H

#include <CGAL/basic.h>

#include <CGAL/leda_bigfloat.h>
#include <CGAL/CORE_BigFloat.h>
#include <CGAL/Gmpfr.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {
    
// Don't define default, results in more convinient compiler messages 
template< class Type > class Float_traits;
// {
// public:        
//   typedef Null_functor    Get_mantissa;
//   typedef Null_functor    Get_exponent;  
//   typedef Null_functor    Mul_by_pow_of_2;
// };
    
#ifdef CGAL_USE_LEDA

// Specialization for leda_bigfloat
template<>
class Float_traits< leda_bigfloat > {
public:
  struct Get_mantissa
    : public std::unary_function< leda_bigfloat, leda_integer > {
    leda_integer operator()( const leda_bigfloat& x ) const {
      //std::cout << x.get_significant() << std::endl;
      return x.get_significant();                
    }
  };
        
  struct Get_exponent
    : public std::unary_function< leda_bigfloat, long > {
    long operator()( const leda_bigfloat& x ) const {
      return x.get_exponent().to_long();                
    }
  };

  struct Mul_by_pow_of_2
    : public std::binary_function< leda_bigfloat, long, leda_bigfloat> {
    leda_bigfloat operator()( const leda_bigfloat& a, long e ) const {
      return leda_bigfloat(a.get_significant(), a.get_exponent()+e);
    }
  };
};

#endif    
    
#ifdef CGAL_USE_CORE

// Specialization for CORE::BigFloat
template<>
class Float_traits< CORE::BigFloat > {
public:
      
  struct Get_mantissa
    : public std::unary_function< CORE::BigFloat, CORE::BigInt > {
    CORE::BigInt operator()( const CORE::BigFloat& x ) const { 
      return x.m();
    }
  };
        
  struct Get_exponent
    : public std::unary_function< CORE::BigFloat, long > {
    long operator()( const CORE::BigFloat& x ) const {
      return 14*x.exp(); // The basis is 8092                 
    }
  };

  struct Mul_by_pow_of_2
    : public std::binary_function
    < CORE::BigFloat, long , CORE::BigFloat> {
    CORE::BigFloat operator()( const CORE::BigFloat& a, long e ) const {
      return a*CORE::BigFloat::exp2(e);
    }
  };

};
#endif    

template<> class Float_traits< Gmpfr > {
public:
  struct Get_mantissa
    : public std::unary_function< Gmpfr, Gmpz > {
    Gmpz operator()( const Gmpfr& x ) const {

      //std::cout << "Get_mantissa" <<std::endl;

      std::pair<Gmpz,long> pair(x.to_integer_exp());
      Gmpfr tmp (pair.first, x.get_prec());
      if (pair.second > 0)
        mpfr_mul_2ui(tmp.fr(),tmp.fr(),pair.second,GMP_RNDN);
      else 
        mpfr_div_2ui(tmp.fr(),tmp.fr(),-pair.second,GMP_RNDN);
      CGAL_postcondition(x == tmp);
      return pair.first;  
    }
  };
  
  struct Get_exponent
    : public std::unary_function< Gmpfr, long > {
    long operator()( const Gmpfr& x ) const {

      //std::cout << "Get_exponent" <<std::endl;

      std::pair<Gmpz,long> pair(x.to_integer_exp());
      Gmpfr tmp (pair.first, x.get_prec());
      if (pair.second > 0)
        mpfr_mul_2ui(tmp.fr(),tmp.fr(),pair.second,GMP_RNDN);
      else 
        mpfr_div_2ui(tmp.fr(),tmp.fr(),-pair.second,GMP_RNDN);
      CGAL_postcondition(x == tmp);
      return pair.second;
    }
  };
    
struct Mul_by_pow_of_2
  : public std::binary_function< Gmpfr, Gmpz, Gmpfr> {
  Gmpfr operator()( const Gmpfr& a, long e ) const {

    //std::cout << "Mul_by_pow_of_2" <<std::endl;

    Gmpfr result; 
    if (e >= 0 ){
      mpfr_mul_2si (result.fr(), a.fr(), (unsigned long)  e, mpfr_get_default_rounding_mode());
      CGAL_postcondition(a * CGAL::ipower(Gmpfr(2),e) == result);
    }
    else{
      mpfr_div_2si (result.fr(), a.fr(), (unsigned long) -e, mpfr_get_default_rounding_mode());
      CGAL_postcondition(a / CGAL::ipower(Gmpfr(2),-e) != result);
    }
    return result;
  }
};
};
} //namespace CGALi



CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H
