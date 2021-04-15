// Copyright (c) 1999,2001,2007
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
      : public CGAL::cpp98::binary_function< Type, Type&,
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
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          return Interval_nt<true>(x).pair();
        }
    };
};

// unsigned long long
template <> class Real_embeddable_traits< unsigned long long >
  : public INTERN_RET::Real_embeddable_traits_base< unsigned long long , CGAL::Tag_true > {
  public:

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          return Interval_nt<true>(x).pair();
        }
    };
};

#ifdef BOOST_HAS_INT128
// __int128
template<> class Algebraic_structure_traits< boost::int128_type >
  : public Algebraic_structure_traits_base< boost::int128_type,
                                            Euclidean_ring_tag > {

  public:
    typedef Tag_true            Is_exact;
    typedef Tag_false           Is_numerical_sensitive;

    typedef INTERN_AST::Div_per_operator< Type >  Div;
    typedef INTERN_AST::Mod_per_operator< Type >  Mod;

};

template <> class Real_embeddable_traits< boost::int128_type >
  : public INTERN_RET::Real_embeddable_traits_base< boost::int128_type , CGAL::Tag_true > {
  public:

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          return (Interval_nt<>((double)x)+Interval_nt<>::smallest()).pair();
        }
    };
};

// unsigned __int128
template <> class Real_embeddable_traits< boost::uint128_type >
  : public INTERN_RET::Real_embeddable_traits_base< boost::uint128_type , CGAL::Tag_true > {
  public:

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
          return (Interval_nt<>((double)x)+Interval_nt<>::smallest()).pair();
        }
    };
};
#endif


} //namespace CGAL

#include <CGAL/Interval_nt.h>

#endif // CGAL_LONG_LONG_H
