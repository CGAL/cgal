// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_REAL_EMBEDDABLE_TRAITS_H
#define CGAL_SQRT_EXTENSION_REAL_EMBEDDABLE_TRAITS_H

#include <CGAL/basic.h>

namespace CGAL {

template< class COEFF, class ROOT, class ACDE_TAG, class FP_TAG >
class Real_embeddable_traits< Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG> >
  : public INTERN_RET::Real_embeddable_traits_base<
                  Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG>,
                  typename Real_embeddable_traits<COEFF>::Is_real_embeddable > {
  public:
  typedef Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG> Type;

    class Sgn
        : public CGAL::cpp98::unary_function< Type, ::CGAL::Sign >{
    public:
        ::CGAL::Sign operator()( const Type& x ) const {
            return x.sign();
        }
    };

    class Compare
        : public CGAL::cpp98::binary_function< Type, Type, Comparison_result > {
    public:
        Comparison_result operator()( const Type& x, const Type& y) const {
            // must be from the same extension
            return x.compare(y);
        }
        Comparison_result operator()( const Type& x, const COEFF& y) const {
            // must be from the same extension
            return x.compare(y);
        }
        Comparison_result operator()( const COEFF& x, const Type& y) const {
            // must be from the same extension
            return CGAL::opposite( y.compare(x) );
        }

        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT( Type,
                Comparison_result )
    };

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double,double> operator()(const Type& x) const {
            return x.to_interval();
        }
    };

    class To_double
      : public CGAL::cpp98::unary_function< Type, double > {
      public:
        // The main problem here is, that even tough the total
        // expression fits into double, one of the coefficients
        // or the root may not. ?? !
        double operator()(const Type& x) const {
            if(x.is_extended()){
                return to_double(x.a0()) +  to_double(x.a1())
                    * (std::sqrt) (to_double(x.root()));
            }else{
                return CGAL_NTS to_double(x.a0());
            }
        }
    };
};
} //namespace CGAL

#endif
