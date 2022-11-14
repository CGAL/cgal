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

#ifndef CGAL_GMPZF_H
#define CGAL_GMPZF_H

// includes
#include <CGAL/number_type_basic.h>
#include <CGAL/Gmp_coercion_traits.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Interval_nt.h>

namespace CGAL {

// Algebraic structure traits
template <> class Algebraic_structure_traits< Gmpzf >
    : public Algebraic_structure_traits_base< Gmpzf, Euclidean_ring_tag >  {
public:
    typedef Tag_true            Is_exact;

    struct Is_zero
        : public CGAL::cpp98::unary_function< Type, bool > {
    public:
        bool operator()( const Type& x ) const {
            return x.is_zero();
        }
    };

    struct Integral_division
        : public CGAL::cpp98::binary_function< Type,
                                Type,
                                Type > {
    public:
        Type operator()(
                const Type& x,
                const Type& y ) const {
            return x.integral_division(y);
        }
    };

    struct Gcd
        : public CGAL::cpp98::binary_function< Type,
                                Type,
                                Type > {
    public:
        Type operator()(
                const Type& x,
                const Type& y ) const {
            return x.gcd(y);
        }
        CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR(int)
    };

    class Div
        : public CGAL::cpp98::binary_function< Type, Type, Type > {
    public:
        Type operator()( const Type& x, const Type& y ) const {
            return Type(x).div( y );
        }
    };

    typedef INTERN_AST::Mod_per_operator< Type > Mod;

  class Is_square
    : public CGAL::cpp98::binary_function< Type, Type&, bool > {
  public:
    bool operator()( const Type& x, Type& y ) const {
      y = CGAL::approximate_sqrt(x);
      return y * y == x;
    }
    bool operator()( const Type& x) const {
      Type dummy;
      return operator()(x,dummy);
    }
  };
};


// Real embeddable traits
template <>
class Real_embeddable_traits< Gmpzf >
    : public INTERN_RET::Real_embeddable_traits_base< Gmpzf , CGAL::Tag_true > {

    typedef Algebraic_structure_traits<Gmpzf> AST;
public:
  typedef AST::Is_zero Is_zero;

    struct Sgn
        : public CGAL::cpp98::unary_function< Type, ::CGAL::Sign > {
    public:
        ::CGAL::Sign operator()( const Type& x ) const {
            return x.sign();
        }
    };

    struct Compare
        : public CGAL::cpp98::binary_function< Type,
                                  Type,
                                  Comparison_result > {
    public:
        Comparison_result operator()(
                const Type& x,
                const Type& y ) const {
            return x.compare(y);
        }
    };

    struct To_double
        : public CGAL::cpp98::unary_function< Type, double > {
    public:
        double operator()( const Type& x ) const {
            return x.to_double();
        }
    };

    struct To_interval
        : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
    public:
        std::pair<double, double> operator()( const Type& x ) const {
            return x.to_interval();
        }
    };
};

// specialization of to double functor
template<>
class Real_embeddable_traits< Quotient<Gmpzf> >
    : public
INTERN_QUOTIENT::Real_embeddable_traits_quotient_base< Quotient<Gmpzf> >
{
public:
    struct To_double: public CGAL::cpp98::unary_function<Quotient<Gmpzf>, double>{
        inline
        double operator()(const Quotient<Gmpzf>& q) const {
          std::pair<double, long> n = q.numerator().to_double_exp();
          std::pair<double, long> d = q.denominator().to_double_exp();
          double scale = std::ldexp(1.0,
                                    static_cast<int>(n.second - d.second));
          return (n.first / d.first) * scale;
        }
    };
    struct To_interval
        : public CGAL::cpp98::unary_function<Quotient<Gmpzf>, std::pair<double,double> >{
        inline
        std::pair<double,double> operator()(const Quotient<Gmpzf>& q) const {
          // do here as MP_Float does
          std::pair<std::pair<double, double>, long> n =
            q.numerator().to_interval_exp();
          std::pair<std::pair<double, double>, long> d =
            q.denominator().to_interval_exp();

          CGAL_assertion_msg(CGAL::abs(double(n.second) - double(d.second)) < (1<<30)*2.0,
                     "Exponent overflow in Quotient<MP_Float> to_interval");
          return ldexp(Interval_nt<>(n.first) / Interval_nt<>(d.first),
                       static_cast<int>(n.second - d.second)).pair();
        }
    };
};

} //namespace CGAL

//since types are included by Gmp_coercion_traits.h:
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>

#endif // CGAL_GMPZF_H

// ===== EOF ==================================================================
