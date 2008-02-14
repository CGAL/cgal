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

#ifndef CGAL_ALGEBRAIC_KERNEL_D_REAL_EMBEDDABLE_EXTENSION_H
#define CGAL_ALGEBRAIC_KERNEL_D_REAL_EMBEDDABLE_EXTENSION_H

#include <CGAL/basic.h>

#ifdef CGAL_USE_LEDA

#include <LEDA/numbers/digit.h>
// #include <CGAL/Algebraic_kernel_d/leda_interval_support.h>
#include <CGAL/leda_interval_support.h>

#endif

#ifdef CGAL_USE_CORE

// #include <CGAL/Number_types/core_interval_support.h>
#include <CGAL/core_interval_support.h>

#endif

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// TODO: Implement array in source code file
//    extern const signed char floor_log2_4bit[16]; // see src/floor_log2_4bit.C
    
    template< class Type >
    class Real_embeddable_extension {
      public:
        typedef Null_functor Ceil_log2_abs;
        typedef Null_functor Floor_log2_abs;
    };
    
    // Functor adapting functions
    template< class NT >
    long floor_log2_abs( NT x ) {
        return typename Real_embeddable_extension< NT >::Floor_log2_abs()( x );
    }
    
    template< class NT >
    long ceil_log2_abs( NT x ) {
        return typename Real_embeddable_extension< NT >::Ceil_log2_abs()( x );
    }

    // Specialization for long
    template<>
    class Real_embeddable_extension< long > {            
      public:      
        struct Ceil_log2_abs
            : public Unary_function< long, long > {
            long operator()( long x ) {
                if (x < 0) x = -x;
                CGAL_precondition(x > 0);
                if (x == 1) return 0;
                return Floor_log2_abs()(x-1) + 1;
            }
        };

        struct Floor_log2_abs
            : public Unary_function< long, long > {
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
    };


#ifdef CGAL_USE_LEDA
    // Specialization for leda_integer
    
    template<>
    class Real_embeddable_extension< leda_integer > {
      public:
        struct Ceil_log2_abs
            : public Unary_function< leda_integer, long > {
                long operator()( leda_integer x ) const {
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
    };
    
    template<>
    class Real_embeddable_extension< leda_bigfloat > {
      public:
        struct Floor_log2_abs
            : public Unary_function< leda_bigfloat, long > {
            long operator()( leda_bigfloat x ) const {
                CGAL_precondition(CGAL::sign(x) != CGAL::ZERO);
                ::leda::integer abs_sign = abs(x.get_significant());
                return (x.get_exponent() + ::leda::log(abs_sign)).to_long();
                
            }
        };
        
        struct Ceil_log2_abs
            : public Unary_function< leda_bigfloat, long > {
            long operator()( leda_bigfloat x ) const {
                CGAL_precondition(CGAL::sign(x) != CGAL::ZERO);
                return ::leda::ilog2(x).to_long();                
            }
        };
    };
    
    template<>
    class Real_embeddable_extension< leda_bigfloat_interval > {
        public:
            struct Floor_log2_abs
                : public Unary_function< leda_bigfloat_interval, long > {
                  
                  result_type operator() (argument_type x) const {
                    CGAL_precondition(! ::boost::numeric::in_zero(x));
                    return CGALi::floor_log2_abs(::boost::numeric::abs(x).lower());
                  }                    
             };
             
            struct Ceil_log2_abs
                : public Unary_function< leda_bigfloat_interval, long > {
                long operator()( leda_bigfloat_interval x ) const {
                    CGAL_precondition(!(::boost::numeric::in_zero(x) && 
                                        ::boost::numeric::singleton(x)));
                    return CGALi::ceil_log2_abs(::boost::numeric::abs(x).upper());                    
                }
            };
    };
                  
#endif
    
#ifdef CGAL_USE_CORE

    // Specialization for CORE::BigInt
    template<>
    class Real_embeddable_extension< CORE::BigInt > {
      public:
        struct Floor_log2_abs
            : public Unary_function< CORE::BigInt, long > {
            long operator()( CORE::BigInt x ) const {
                return CORE::floorLg(x);
            }            
        };
        
        struct Ceil_log2_abs
            : public Unary_function< CORE::BigInt, long > {
            long operator()( CORE::BigInt x ) const {
                return CORE::ceilLg(x);
            }
        };
    };

    // Specialization for CORE::BigFloat
    template<>
    class Real_embeddable_extension< CORE::BigFloat > {
      public:
        struct Floor_log2_abs
            : public Unary_function< CORE::BigFloat, long > {
            long operator()( CORE::BigFloat x ) const {
                CGAL_precondition(!CGAL::in_zero(x));
                x = CGAL::abs(x);
                return CORE::floorLg(x.m()-x.err())+x.exp()*14;
            }            
        };
        
        struct Ceil_log2_abs
            : public Unary_function< CORE::BigFloat, long > {
            long operator()( CORE::BigFloat x ) const {
                // (already commented out in EXACUS)...
                //   NiX_precond(!(NiX::in_zero(x) && NiX::singleton(x)));
                x = CGAL::abs(x);
                return CORE::ceilLg(x.m()+x.err())+x.exp()*14;
            }
        };
    };

#endif
        
} //namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_REAL_EMBEDDABLE_EXTENSION_H
