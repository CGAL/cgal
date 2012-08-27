// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================

// TODO: Some comments are original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_ALGEBRAIC_KERNEL_D_REAL_EMBEDDABLE_EXTENSION_H
#define CGAL_ALGEBRAIC_KERNEL_D_REAL_EMBEDDABLE_EXTENSION_H

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_d/Float_traits.h>

#ifdef CGAL_USE_LEDA

#include <CGAL/leda_integer.h>
#include <CGAL/leda_bigfloat.h>
#include <CGAL/leda_bigfloat_interval.h>
#include <LEDA/numbers/digit.h>
#endif

#ifdef CGAL_USE_CORE
#include <CGAL/CORE_BigFloat.h>
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#endif 

#ifdef CGAL_USE_MPFR
#include <CGAL/Gmpfr.h>
#endif 

#ifdef CGAL_USE_MPFI
#include <CGAL/Gmpfi.h>
#endif 


namespace CGAL {

namespace internal {

// TODO: Implement array in source code file
//    extern const signed char floor_log2_4bit[16]; // see src/floor_log2_4bit.C
    
// Don't define default, results in more convinient compiler messages 
template< class Type > class Real_embeddable_extension;
// {
// public:
//   typedef Null_functor Ceil_log2_abs;
//   typedef Null_functor Floor_log2_abs;
//   typedef Null_functor Floor;
//   typedef Null_functor Ceil;
// };


    
// Functor adapting functions
template< class NT >
long floor_log2_abs( const NT& x ) {
  return typename Real_embeddable_extension< NT >::Floor_log2_abs()( x );
}
    
template< class NT >
long ceil_log2_abs( const NT& x ) {
  return typename Real_embeddable_extension< NT >::Ceil_log2_abs()( x );
}

template< class NT >
typename Real_embeddable_extension<NT>::Floor::result_type
floor (const NT& x) {
  return typename Real_embeddable_extension<NT>::Floor() (x);
}

template< class NT >
typename Real_embeddable_extension<NT>::Ceil::result_type
ceil (const NT& x) {
  return typename Real_embeddable_extension<NT>::Ceil() (x);
}

// Specialization for long
template<>
class Real_embeddable_extension< long > {            
public:      
  struct Ceil_log2_abs
    : public std::unary_function< long, long > {
    long operator()( long x ) {
      if (x < 0) x = -x;
      CGAL_precondition(x > 0);
      if (x == 1) return 0;
      return Floor_log2_abs()(x-1) + 1;
    }
  };

  struct Floor_log2_abs
    : public std::unary_function< long, long > {
  private:          
    signed char floor_log2_4bit[16];
  public:
    Floor_log2_abs() {
      floor_log2_4bit[ 0] = -42;
      floor_log2_4bit[ 1] = 0;
      floor_log2_4bit[ 2] = 1;
      floor_log2_4bit[ 3] = 1;
      floor_log2_4bit[ 4] = 2;
      floor_log2_4bit[ 5] = 2;
      floor_log2_4bit[ 6] = 2;
      floor_log2_4bit[ 7] = 2;
      floor_log2_4bit[ 8] = 3;
      floor_log2_4bit[ 9] = 3;
      floor_log2_4bit[10] = 3;
      floor_log2_4bit[11] = 3;
      floor_log2_4bit[12] = 3;
      floor_log2_4bit[13] = 3;
      floor_log2_4bit[14] = 3;
      floor_log2_4bit[15] = 3;                            
    }
            
    long operator()( long x ) {
      if (x < 0) x = -x;
      CGAL_precondition(x > 0);
      result_type l = 0;
      while (x > 0xFFFF) { l += 16; x >>= 16; }
      if    (x > 0xFF)   { l +=  8; x >>=  8; }
      if    (x > 0xF)    { l +=  4; x >>=  4; }
      CGAL_assertion(x > 0 && x < 16);
      return l + int(floor_log2_4bit[x]);
    }
  };

  struct Floor
    : public std::unary_function< long, long > {
    long operator() (long x) { return x;}
  };
  struct Ceil
    : public std::unary_function< long, long > {
    long operator() (long x) { return x;}
  };
};


#ifdef CGAL_USE_LEDA
// Specialization for leda_integer
    
template<>
class Real_embeddable_extension< leda_integer > {
public:
  typedef leda_integer Type;

  struct Ceil_log2_abs
    : public std::unary_function< leda_integer, long > {
    long operator()( const leda_integer& x ) const {
      CGAL_precondition(x != leda_integer(0));
      ::leda::digit_sz ldgzeros = ::leda::digLeadingZeros(x.highword());
      result_type l =
        x.used_words() * ::leda::DIGIT_LENGTH - 1 - ldgzeros;
      // look if additional 1-bits force to round up
      ::leda::digit h = 1;
      h <<= ::leda::DIGIT_LENGTH - 1 - ldgzeros;
      int i = x.used_words() - 1;
      CGAL_assertion(x.contents(i) >= h);
      if (x.contents(i) > h) return l+1;
      while (--i >= 0) {
        if (x.contents(i) != 0) return l+1;
      }
      return l;                    
    }
  };

  struct Floor_log2_abs
    : public std::unary_function< leda_integer, long > {
    long operator()( const leda_integer& x ) const {
      CGAL_precondition(x != leda_integer(0));
      ::leda::digit_sz ldgzeros 
          = ::leda::digLeadingZeros(x.highword());
      result_type l =
        x.used_words() * ::leda::DIGIT_LENGTH - 1 - ldgzeros;
      return l;
    }
  };

  struct Floor
    : public std::unary_function< leda_integer, leda_integer > {
    leda_integer operator() (const leda_integer& x) const { return x;}
  };
  struct Ceil
    : public std::unary_function< leda_integer, leda_integer > {
    leda_integer operator() (const leda_integer& x) const { return x;}
  };
};
    
template<>
class Real_embeddable_extension< leda_bigfloat > {
public:

  typedef leda_bigfloat Type;

  struct Floor_log2_abs
    : public std::unary_function< leda_bigfloat, long > {
    long operator()( const leda_bigfloat& x ) const {
      CGAL_precondition(CGAL::sign(x) != CGAL::ZERO);
      ::leda::integer abs_sign = abs(x.get_significant());
      return (x.get_exponent() + ::leda::log(abs_sign)).to_long();
                
    }
  };
        
  struct Ceil_log2_abs
    : public std::unary_function< leda_bigfloat, long > {
    long operator()( const leda_bigfloat& x ) const {
      CGAL_precondition(CGAL::sign(x) != CGAL::ZERO);
      return ::leda::ilog2(x).to_long();                
    }
  };

  struct Floor
    : public std::unary_function< leda_bigfloat, leda_integer > {
    leda_integer operator() ( const leda_bigfloat& x ) const { 
      return leda::to_integer( x, leda::TO_N_INF );
    }
  };
        
  struct Ceil
    : public std::unary_function< leda_bigfloat, leda_integer > {
    leda_integer operator() ( const leda_bigfloat& x ) const { 
      return leda::to_integer( x, leda::TO_P_INF );
    }
  };

};
    
template<>
class Real_embeddable_extension< leda_bigfloat_interval > {
public:
  typedef leda_bigfloat_interval Type;

  struct Floor_log2_abs
    : public std::unary_function< leda_bigfloat_interval, long > {
              
    result_type operator() (const argument_type& x) const {
      CGAL_precondition(! ::boost::numeric::in_zero(x));
      return internal::floor_log2_abs(::boost::numeric::abs(x).lower());
    }                    
  };
        
  struct Ceil_log2_abs
    : public std::unary_function< leda_bigfloat_interval, long > {
    long operator()( const leda_bigfloat_interval& x ) const {
      CGAL_precondition(!(::boost::numeric::in_zero(x) && 
              ::boost::numeric::singleton(x)));
      return internal::ceil_log2_abs(::boost::numeric::abs(x).upper());                    
    }
  };

  struct Floor
    : public std::unary_function< leda_bigfloat_interval, leda_integer > {
    leda_integer operator() ( const leda_bigfloat_interval& x ) 
      const { 
      return internal::floor( x.lower() );
    }
  };
        
  struct Ceil
    : public std::unary_function< leda_bigfloat_interval, leda_integer > {
    leda_integer operator() ( const leda_bigfloat_interval& x ) 
      const { 
      return internal::ceil( x.upper() );
    }
  };
};
                  
#endif
    
#ifdef CGAL_USE_CORE

// Specialization for CORE::BigInt
template<>
class Real_embeddable_extension< CORE::BigInt > {
public:
  typedef CORE::BigInt Type;
  struct Floor_log2_abs
    : public std::unary_function< CORE::BigInt, long > {
    long operator()( const CORE::BigInt& x ) const {
      return CORE::floorLg(x);
    }            
  };
        
  struct Ceil_log2_abs
    : public std::unary_function< CORE::BigInt, long > {
    long operator()( const CORE::BigInt& x ) const {
      return CORE::ceilLg(x);
    }
  };

  struct Floor
    : public std::unary_function< CORE::BigInt, CORE::BigInt > {
    CORE::BigInt operator() (const CORE::BigInt& x) const { 
      return x;
    }
  };
  struct Ceil
    : public std::unary_function< CORE::BigInt, CORE::BigInt > {
    CORE::BigInt operator() (const CORE::BigInt& x) const { 
      return x;
    }
  };
};

// Specialization for CORE::BigFloat
template<>
class Real_embeddable_extension< CORE::BigFloat > {
public:
  typedef CORE::BigFloat Type;
  struct Floor_log2_abs
    : public std::unary_function< CORE::BigFloat, long > {
    long operator()( CORE::BigFloat x ) const {
      CGAL_precondition(!CGAL::zero_in(x));
      x = CGAL::abs(x);
      return CORE::floorLg(x.m()-x.err())+x.exp()*CORE::CHUNK_BIT;
    }            
  };
        
  struct Ceil_log2_abs
    : public std::unary_function< CORE::BigFloat, long > {
    long operator()( CORE::BigFloat x ) const {
      // (already commented out in EXACUS)...
      //   NiX_precond(!(NiX::in_zero(x) && NiX::singleton(x)));
      x = CGAL::abs(x);
      return CORE::ceilLg(x.m()+x.err())+x.exp()*CORE::CHUNK_BIT;
    }
  };

  struct Floor
    : public std::unary_function< CORE::BigFloat, CORE::BigInt > {
    CORE::BigInt operator() ( const CORE::BigFloat& x ) const { 
      CORE::BigInt xi = x.BigIntValue();
      if(x.sign() < 0 && x.cmp(xi)!=0) {
        xi--;
      }
      return xi;
    }
  };
        
  struct Ceil
    : public std::unary_function< CORE::BigFloat, CORE::BigInt > {
    CORE::BigInt operator() ( const CORE::BigFloat& x ) const { 
      CORE::BigInt xi = x.BigIntValue();
      if(x.sign() >0 && x.cmp(xi)!=0) {
        xi++;
      }
      return xi;
    }
  };

};

#endif // CORE

#if CGAL_USE_GMP

// Specialization for Gmpz
template<>
class Real_embeddable_extension< Gmpz > {
public:
  typedef Gmpz Type;

  struct Floor_log2_abs
    : public std::unary_function< Gmpz, long > {
    long operator()( const Gmpz& x ) const {
      CGAL_precondition(!CGAL::is_zero(x));
      return mpz_sizeinbase(x.mpz(),2)-1;
    }            
  };
        
  struct Ceil_log2_abs
    : public std::unary_function< Gmpz, long > {
    long operator()( const Gmpz& x ) const {
      long pos  = mpz_scan1(x.mpz(),0);
      long size = mpz_sizeinbase(x.mpz(),2);
      if (pos == size-1) 
        return size-1;
      else 
        return size;
    }
  };

  struct Floor
    : public std::unary_function< Gmpz, Gmpz > {
    Gmpz operator() (const Gmpz& x) const { 
      return x;
    }
  };
  struct Ceil
    : public std::unary_function< Gmpz, Gmpz > {
    Gmpz operator() (const Gmpz& x) const { 
      return x;
    }
  };
};

#endif 

#ifdef CGAL_USE_MPFR 
template<>
class Real_embeddable_extension< Gmpfr > {
public:
  typedef Gmpfr Type;

  struct Floor_log2_abs
    : public std::unary_function< Gmpfr, long > {
    long operator()( const Gmpfr& x ) const {
      Float_traits<Gmpfr>::Get_mantissa get_mantissa; 
      Float_traits<Gmpfr>::Get_exponent get_exponent; 
      CGAL_precondition(!CGAL::is_zero(x));
      Real_embeddable_extension<Gmpz>::Floor_log2_abs floor_log2_abs;
      return floor_log2_abs(get_mantissa(x))+get_exponent(x);
    }
  };
        
  struct Ceil_log2_abs
    : public std::unary_function< Gmpfr, long > {
    long operator()( const Gmpfr& x ) const {
      Float_traits<Gmpfr>::Get_mantissa get_mantissa; 
      Float_traits<Gmpfr>::Get_exponent get_exponent; 
      CGAL_precondition(!CGAL::is_zero(x));
      Real_embeddable_extension<Gmpz>::Ceil_log2_abs ceil_log2_abs;
      return ceil_log2_abs(get_mantissa(x))+get_exponent(x);
    }
  };

  struct Floor
    : public std::unary_function< Gmpfr, Gmpz > {
    Gmpz operator() ( const Gmpfr& x ) const {  
      Gmpz result; 
      mpfr_get_z (result.mpz(),x.fr(),GMP_RNDD);
      return result; 
    }
  };
        
  struct Ceil
    : public std::unary_function< Gmpfr, Gmpz > {
    Gmpz operator() ( const Gmpfr& x ) const { 
      Gmpz result; 
      mpfr_get_z (result.mpz(),x.fr(),GMP_RNDU);
      return result; 
    }
  };

};
#endif     

#ifdef CGAL_USE_MPFI
template<>
class Real_embeddable_extension< Gmpfi > {
public:
  typedef Gmpfi Type;

  struct Floor_log2_abs
    : public std::unary_function< Gmpfi, long > {
    result_type operator() (const argument_type& x) const {
      CGAL_precondition(!x.is_zero());
      return internal::floor_log2_abs(x.abs().inf());
    }                    
  };
        
  struct Ceil_log2_abs
    : public std::unary_function< Gmpfi, long > {
    long operator()( const Gmpfi& x ) const {
      CGAL_precondition(!x.inf().is_zero() || !x.sup().is_zero());
      return internal::ceil_log2_abs(x.abs().sup());                    
    }
  };

  struct Floor
    : public std::unary_function< Gmpfi, Gmpz > {
    Gmpz operator() ( const Gmpfi& x ) 
      const { 
      return internal::floor( x.inf() );
    }
  };
        
  struct Ceil
    : public std::unary_function< Gmpfi, Gmpz > {
    Gmpz operator() ( const Gmpfi& x ) 
      const { 
      return internal::ceil( x.sup() );
    }
  };
};
  
#endif
        
} //namespace internal

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_D_REAL_EMBEDDABLE_EXTENSION_H
