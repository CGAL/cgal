// Copyright (c) 1999,2001,2007  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Stefan Schirra, Michael Hemmer

// ISO C++ does not support `long long', but ISO C does, which means the next
// revision of ISO C++ probably will too.  However, currently, g++ -pedantic
// produces a warning so we don't include this file by default.

#ifndef CGAL_LONG_LONG_H
#define CGAL_LONG_LONG_H

#include <CGAL/number_type_basic.h>

namespace CGAL {

template<> class Algebraic_structure_traits< long long int >
  : public Algebraic_structure_traits_base< long long int,
                                            Euclidean_ring_tag > {

  public:
    typedef Tag_true            Is_exact;
    typedef Tag_false           Is_numerical_sensitive;

    typedef INTERN_AST::Div_per_operator< Type >  Div;
    typedef INTERN_AST::Mod_per_operator< Type >  Mod;

    class Is_square
      : public std::binary_function< Type, Type&,
                                bool > {
      public:
        bool operator()( const Type& x,
                         Type& y ) const {
          y = (Type) std::sqrt( (double)x );
          return x == y * y;
        }
        bool operator()( const Type& x) const {
            Type y
                = (Type) std::sqrt( (double)x );
          return x == y * y;
        }
    };
};

template <> class Real_embeddable_traits< long long int >
  : public INTERN_RET::Real_embeddable_traits_base< long long int , CGAL::Tag_true > {
  public:

    class To_interval
      : public std::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
          Interval_nt<false> approx ((double) x);
          FPU_set_cw(CGAL_FE_UPWARD);
          approx += Interval_nt<false>::smallest();
          return approx.pair();
        }
    };
};

} //namespace CGAL

#include <CGAL/Interval_nt.h>

#endif // CGAL_LONG_LONG_H
