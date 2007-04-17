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

#ifndef CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H

#include <CGAL/basic.h>


CGAL_BEGIN_NAMESPACE

namespace CGALi {
    
    template< class Type >
    class Float_traits {
      public:
        
        typedef Null_functor    Get_mantissa;
        typedef Null_functor    Get_exponent;  
    };
    
#ifdef CGAL_USE_LEDA

    // Specialization for leda_bigfloat
    template<>
    class Float_traits< leda_bigfloat > {
      public:
      
        struct Get_mantissa
            : public Unary_function< leda_bigfloat, leda_integer > {
            leda_integer operator()( const leda_bigfloat& x ) const {
                return x.get_significant();                
            }
        };
        
        struct Get_exponent
            : public Unary_function< leda_bigfloat, long > {
            long operator()( const leda_bigfloat& x ) const {
                return x.get_exponent().to_long();                
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
            : public Unary_function< CORE::BigFloat, CORE::BigInt > {
            CORE::BigInt operator()( const CORE::BigFloat& x ) const {
                return x.m();
            }
        };
        
        struct Get_exponent
            : public Unary_function< CORE::BigFloat, long > {
            long operator()( const CORE::BigFloat& x ) const {
                return 14*x.exp(); // The basis is 8092                 
            }
        };
    };

#endif    
    
} //namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_FLOAT_TRAITS_H
