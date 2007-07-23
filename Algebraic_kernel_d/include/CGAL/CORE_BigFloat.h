
// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Number_types/include/CGAL/CORE_BigFloat.h $
// $Id: CORE_BigFloat.h 38140 2007-04-16 08:57:45Z hemmer $
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>

#ifndef CGAL_CORE_BIGFLOAT_H
#define CGAL_CORE_BIGFLOAT_H

#include <CGAL/number_type_basic.h>
#include <CGAL/CORE_coercion_traits.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi{
// Forward declarations of core_interval_support.h
long inline get_significant_bits(const CORE::BigFloat& x);
CORE::BigFloat inline round(const CORE::BigFloat& x, long rel_prec);
long inline get_significant_bits(const CORE::BigFloat& x);
long inline set_precision(CORE::BigFloat, long prec);
long inline get_precision(CORE::BigFloat);
CORE::BigFloat inline upper(CORE::BigFloat x);
CORE::BigFloat inline lower(CORE::BigFloat x);
bool inline in_zero(CORE::BigFloat x);
bool inline overlap(CORE::BigFloat x, CORE::BigFloat y);
CORE::BigFloat inline hull(CORE::BigFloat x, CORE::BigFloat y);
bool inline singleton(CORE::BigFloat x);
CORE::BigFloat inline convert_to_bfi(const ::CORE::Expr& x);
CORE::BigFloat inline convert_to_bfi(const ::CORE::BigInt& x);
CORE::BigFloat inline convert_to_bfi(const ::CORE::BigRat& x);
}

//
// Algebraic structure traits
//
template <> class Algebraic_structure_traits< CORE::BigFloat >
  : public Algebraic_structure_traits_base< CORE::BigFloat,
                                            Field_with_kth_root_tag >  {
    typedef CORE::BigFloat NT;
  public:
    typedef Tag_false          Is_exact;
    typedef Tag_true           Is_numerical_sensitive;

    class Sqrt
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& a ) const {
            // What I want is a sqrt computed with ::CORE::defRelPrec bits.
            // And not ::CORE::defBFsqrtAbsPrec as CORE does. 
            
            CGAL_precondition(::CORE::defRelPrec.toLong() > 0);
            CGAL_precondition(a > 0);
            
            NT x = CGAL::CGALi::round(a, ::CORE::defRelPrec.toLong()*2);
            CGAL_postcondition(x > 0); 

            NT tmp1 = 
                CORE::BigFloat(x.m(),0,0).sqrt(::CORE::defRelPrec.toLong());
            NT err  =  
                NT(0,long(std::sqrt(double(x.err()))),0) 
                * CORE::BigFloat::exp2(x.exp()*7);
            NT result = tmp1*CORE::BigFloat::exp2(x.exp()*7) + err;
           
            CGAL_postcondition(result >= 0);
#ifndef NDEBUG
            NT tmp = result * result; 
            //bfi_traits::Upper upper;
            //bfi_traits::Lower lower;
            CGAL_postcondition(CGAL::CGALi::lower(tmp) <= CGAL::CGALi::lower(a));
            CGAL_postcondition(CGAL::CGALi::upper(tmp) >= CGAL::CGALi::upper(a));
#endif
            return result;
        }
    };

    class Kth_root
      : public Binary_function<int, Type, Type> {
      public:
        Type operator()( int k,
                                        const Type& x) const {
            CGAL_precondition_msg( k > 0, "'k' must be positive for k-th roots");
            // CORE::radical isn't implemented for negative values of x, so we
            //  have to handle this case separately
            if( x < 0 && k%2 != 0) {
              return Type(-CORE::radical( -x, k ) );
            }
            return Type( CORE::radical( x, k ) );
        }
    };
};

//
// Real embeddable traits
//
template <> class Real_embeddable_traits< CORE::BigFloat >
  : public Real_embeddable_traits_base< CORE::BigFloat > {
    typedef CORE::BigFloat NT;
  public:

    class Abs
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            return CORE::abs( x );
        }
    };

    class Sign
      : public Unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
            return (::CGAL::Sign) CORE::sign( x );
        }
    };

    class Compare
      : public Binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& x,
                                            const Type& y ) const {
          return (Comparison_result) CORE::cmp( x, y );
        }
    };

    class To_double
      : public Unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          // this call is required to get reasonable values for the double
          // approximation
          return x.doubleValue();
        }
    };

    class To_interval
      : public Unary_function< Type, std::pair< double, double > > {
    public:
        std::pair<double, double> operator()( const Type& x ) const {
            double lb,ub;
           
            NT x_lower = CGAL::CGALi::lower(CGAL::CGALi::round(CGAL::CGALi::lower(x),52));
            NT x_upper = CGAL::CGALi::upper(CGAL::CGALi::round(CGAL::CGALi::upper(x),52));
            
            // since matissa has 52 bits only, conversion to double is exact 
            lb = x_lower.doubleValue();
            CGAL_postcondition(lb == x_lower);
            ub = x_upper.doubleValue();
            CGAL_postcondition(ub == x_upper);             
            
            result_type result(lb,ub);
            CGAL_postcondition( result.first  <=  CORE::Expr(CGAL::CGALi::lower(x)));
            CGAL_postcondition( result.second >=  CORE::Expr(CGAL::CGALi::upper(x)));
            return result;       
           
        }
    };
};

CGAL_END_NAMESPACE

//since types are included by CORE_coercion_traits.h:
#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_BigFloat.h>

#include <CGAL/Number_types/core_interval_support.h>

#endif // CGAL_CORE_BIGFLOAT_H
