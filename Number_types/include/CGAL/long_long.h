// Copyright (c) 1999,2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
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

#include <CGAL/basic.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/utils.h>
#include <CGAL/functional_base.h> // Unary_function, Binary_function

CGAL_BEGIN_NAMESPACE

template <> struct Number_type_traits<long long int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_true   Has_division;
  typedef Tag_false  Has_sqrt;

  typedef Tag_false  Has_exact_ring_operations;
  typedef Tag_false  Has_exact_division;
  typedef Tag_false  Has_exact_sqrt;
};

template<> class Algebraic_structure_traits< long long int >
  : public Algebraic_structure_traits_base< long long int, 
                                            CGAL::Euclidean_ring_tag > {

  public:
    typedef CGAL::Tag_true            Is_exact;
    
    typedef CGAL::INTERN_AST::Div_per_operator< Algebraic_structure >  Div;
    typedef CGAL::INTERN_AST::Mod_per_operator< Algebraic_structure >  Mod;       

    class Is_square 
      : public Binary_function< Algebraic_structure, Algebraic_structure&,
                                bool > {
      public:
        bool operator()( const Algebraic_structure& x,
                         Algebraic_structure& y ) const {
          y = (Algebraic_structure) CGAL_CLIB_STD::sqrt( (double)x );
          return x == y * y;
        }
        bool operator()( const Algebraic_structure& x) const {
            Algebraic_structure y 
                = (Algebraic_structure) CGAL_CLIB_STD::sqrt( (double)x );
          return x == y * y;
        }
    };
};

template <> class Real_embeddable_traits< long long int > 
  : public Real_embeddable_traits_base< long long int > {
  public:
          
    typedef CGAL::INTERN_RET::To_double_by_conversion< Real_embeddable >
                                                                      To_double;

    class To_interval 
      : public Unary_function< Real_embeddable, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Real_embeddable& x ) const {
          Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
          CGAL::Interval_nt<false> approx ((double) x);
          FPU_set_cw(CGAL_FE_UPWARD);
          approx += Interval_nt<false>::smallest();
          return approx.pair();          
        }
    };
};

#if (defined(__sparc__) || defined(__sparc) || defined(sparc)) || \
    (defined(__sgi__)   || defined(__sgi)   || defined(sgi)) || \
    (defined(__i386__)  || defined(__i386)  || defined(i386)) || \
    (defined(__ppc__)   || defined(__ppc)   || defined(ppc)) || \
    (defined(__powerpc__) || defined(__powerpc) || defined(powerpc))
typedef  long long int           Integer64;
typedef  unsigned long long int  UInteger64;
#define CGAL_HAS_INTEGER64
#endif

CGAL_END_NAMESPACE

#include <CGAL/Interval_nt.h>

#endif // CGAL_LONG_LONG_H
