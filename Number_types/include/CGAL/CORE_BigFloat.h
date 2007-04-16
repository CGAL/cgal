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
          return CORE::sqrt( x );
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
                        
            CORE::Expr lower_expr = ::CORE::BigFloat(x.m()-x.err(),0,x.exp());
            CORE::Expr upper_expr = ::CORE::BigFloat(x.m()-x.err(),0,x.exp());

            std::pair<double, double> lower,upper;
            lower_expr.doubleInterval(lower.first, lower.second);
            upper_expr.doubleInterval(upper.first, upper.second);
            
            CGAL_postcondition(lower.first  <= upper.first);
            CGAL_postcondition(lower.second <= upper.second);

            return std::pair<double, double>(lower.first,upper.second);
        }
    };
};

CGAL_END_NAMESPACE

//since types are included by CORE_coercion_traits.h:
#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_BigFloat.h>

#endif // CGAL_CORE_BIGFLOAT_H
