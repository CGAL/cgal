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
// $URL$
// $Id$
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>

#ifndef CGAL_CORE_BIGFLOAT_H
#define CGAL_CORE_BIGFLOAT_H

#include <CGAL/number_type_basic.h>
#include <CGAL/CORE_coercion_traits.h>

CGAL_BEGIN_NAMESPACE

// forward declarations of interval support
namespace CGALi {
CORE::BigFloat 
inline 
round(const CORE::BigFloat& x, long rel_prec );

CORE::BigFloat 
inline
upper(CORE::BigFloat x);

CORE::BigFloat 
inline 
lower(CORE::BigFloat x);

} // namespace CGALi

//
// Algebraic structure traits
//
template <> class Algebraic_structure_traits< CORE::BigFloat >
  : public Algebraic_structure_traits_base< CORE::BigFloat,
                                            Field_with_kth_root_tag >  {
  public:
    typedef Tag_false          Is_exact;
    typedef Tag_true           Is_numerical_sensitive;

    class Sqrt
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            // What I want is a sqrt computed with ::CORE::defRelPrec bits.
            // And not ::CORE::defBFsqrtAbsPrec as CORE does. 
            
            CGAL_precondition(::CORE::defRelPrec.toLong() > 0);
            CGAL_precondition(x > 0);
            
            Type a = CGALi::round(x, ::CORE::defRelPrec.toLong()*2);
            CGAL_postcondition(a > 0); 

            Type tmp1 = 
                CORE::BigFloat(a.m(),0,0).sqrt(::CORE::defRelPrec.toLong());
            Type err  =  
                Type(0,long(std::sqrt(double(a.err()))),0) 
                * CORE::BigFloat::exp2(a.exp()*7);
            Type result = tmp1*CORE::BigFloat::exp2(a.exp()*7) + err;
           
            CGAL_postcondition(result >= 0);
//#ifndef NDEBUG
//            Type tmp = result * result; 
//            typedef Bigfloat_interval_traits<CORE::BigFloat> bfi_traits;
            //bfi_traits::Upper upper;
            //bfi_traits::Lower lower;
//            CGAL_postcondition(NiX::lower(tmp) <= NiX::lower(x));
//            CGAL_postcondition(NiX::upper(tmp) >= NiX::upper(x));
//#endif
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
  public:

    class Abs
      : public Unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
            Type result; 
          
            if(x.isZeroIn()){
                CORE::BigInt m; 
                if(x.m() < 0 ){
                    m = -(x.m()-x.err());
                }else{
                    m =  x.m()+x.err();
                }
                if(m % 2 == 1) m += 1;
                
                Type upper(m,0,x.exp());
                result = CORE::centerize(CORE::BigFloat(0),upper);
                
                CGAL_postcondition(result.m()-result.err() <= 0); 
                if(result.m()-result.err() != 0){
                    result = this->operator()(result);
                }
                CGAL_postcondition(result.m()-result.err() == 0); 
            }else{
                result = CORE::abs(x);
            }
            CGAL_postcondition(result.m()-result.err() >= 0); 
            CGAL_postcondition(Type(result.m()+result.err(),0,result.exp()) 
                         >= Type(x.m()+x.err(),0,x.exp()));       
            return result;
        }
    };

    class Sign
      : public Unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
            ::CGAL::Sign result =  sign( x.sign());
            return result; 
        }
    };

    class Compare
      : public Binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& x,
                                            const Type& y ) const {
          return (Comparison_result) sign( (x-y).sign());
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
           
            Type x_lower = CGALi::lower(CGALi::round(CGALi::lower(x),52));
            Type x_upper = CGALi::upper(CGALi::round(CGALi::upper(x),52));
            
            // since matissa has 52 bits only, conversion to double is exact 
            lb = x_lower.doubleValue();
            CGAL_postcondition(lb == x_lower);
            ub = x_upper.doubleValue();
            CGAL_postcondition(ub == x_upper);             
            
            std::pair<double, double> result(lb,ub);
            CGAL_postcondition( result.first <=  CORE::Expr(CGALi::lower(x)));
            CGAL_postcondition( result.second >=  CORE::Expr(CGALi::upper(x)));
            return result;      
        }
    };
};

CGAL_END_NAMESPACE

#include <CGAL/Number_types/core_interval_support.h>

//since types are included by CORE_coercion_traits.h:
#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_BigFloat.h>

#endif // CGAL_CORE_BIGFLOAT_H
