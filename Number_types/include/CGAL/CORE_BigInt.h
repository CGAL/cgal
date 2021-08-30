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


#ifndef CGAL_CORE_BIGINT_H
#define CGAL_CORE_BIGINT_H

#include <CGAL/disable_warnings.h>

#include <CGAL/config.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/CORE/BigInt.h>
#include <CGAL/CORE/Expr.h>
#include <CGAL/CORE_coercion_traits.h>

#include <CGAL/Residue.h>
#include <CGAL/Modular_traits.h>

namespace CGAL {

//
// Algebraic structure traits
//
template <> class Algebraic_structure_traits< CORE::BigInt >
  : public Algebraic_structure_traits_base< CORE::BigInt,
                                            Euclidean_ring_tag >  {
  public:
    typedef Tag_true            Is_exact;
    typedef Tag_false           Is_numerical_sensitive;

    typedef INTERN_AST::Is_square_per_sqrt< Type >
                                                                 Is_square;

    typedef INTERN_AST::Div_per_operator< Type > Div;
    typedef INTERN_AST::Mod_per_operator< Type > Mod;

    class Sqrt
      : public CGAL::cpp98::unary_function< Type, Type > {
      public:
        //! computes the largest NT not larger than the square root of \a a.
        Type operator()( const Type& x) const {
          Type result;
          mpz_sqrt(result.get_mp(), x.get_mp());
          return result;
        }
    };


    class Gcd
      : public CGAL::cpp98::binary_function< Type, Type,
                                Type > {
      public:
        Type operator()( const Type& x,
                                        const Type& y) const {
          if ( x == Type(0) && y == Type(0) )
              return Type(0);
          Type result;
          mpz_gcd(result.get_mp(), x.get_mp(), y.get_mp());
          return result;
        }
    };
};

//
// Real embeddable traits
//
template <> class Real_embeddable_traits< CORE::BigInt >
  : public INTERN_RET::Real_embeddable_traits_base< CORE::BigInt , CGAL::Tag_true > {

  public:

    class Abs
      : public CGAL::cpp98::unary_function< Type, Type > {
      public:
        Type operator()( const Type& x ) const {
          return CORE::abs( x );
        }
    };

    class Sgn
      : public CGAL::cpp98::unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
          return (::CGAL::Sign) CORE::sign( x );
        }
    };

    class Compare
      : public CGAL::cpp98::binary_function< Type, Type,
                                Comparison_result > {
      public:
        Comparison_result operator()( const Type& x,
                                            const Type& y ) const {
          return CGAL::sign(::CORE::cmp(x,y));
        }
    };

    class To_double
      : public CGAL::cpp98::unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          // this call is required to get reasonable values for the double
          // approximation
          return x.doubleValue();
        }
    };

    class To_interval
      : public CGAL::cpp98::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x_ ) const {
            CORE::Expr x(x_);
            std::pair<double,double> result;
            x.doubleInterval(result.first, result.second);
            CGAL_expensive_assertion(result.first  <= x);
            CGAL_expensive_assertion(result.second >= x);
            return result;
        }
    };
};

/*! \ingroup NiX_Modular_traits_spec
 *  \brief a model of concept ModularTraits,
 *  specialization of NiX::Modular_traits.
 */
template<>
class Modular_traits< ::CORE::BigInt > {
  typedef Residue RES;
 public:
    typedef ::CORE::BigInt NT;
    typedef CGAL::Tag_true Is_modularizable;
    typedef Residue Residue_type;

    struct Modular_image{
        Residue_type operator()(const NT& a){
            NT tmp = a % NT(RES::get_current_prime());
// TODO: reactivate this assertion
// it fails with core_v1.6x_20040329
//            NiX_assert(tmp.isInt());
            int mi(tmp.longValue());
            if (mi < 0) mi += RES::get_current_prime();
            return Residue_type(mi);
        }
    };
    struct Modular_image_representative{
        NT operator()(const Residue_type& x){
            return NT(x.get_value());
        }
    };
};


template<>
struct Needs_parens_as_product<CORE::BigInt>{
    bool operator()(const CORE::BigInt& x){
        return CGAL_NTS is_negative(x);
    }
};

// Benchmark_rep specialization
template<>
class Benchmark_rep< CORE::BigInt > {
    const CORE::BigInt& t;
public:
    //! initialize with a const reference to \a t.
    Benchmark_rep( const CORE::BigInt& tt) : t(tt) {}
    //! perform the output, calls \c operator\<\< by default.
    std::ostream& operator()( std::ostream& out) const {
            out << t;
            return out;
    }

    static std::string get_benchmark_name() {
        return "Integer";
    }
};


} //namespace CGAL

//since types are included by CORE_coercion_traits.h:
#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_BigRat.h>
#include <CGAL/CORE_BigFloat.h>
#include <CGAL/CORE_arithmetic_kernel.h>

namespace Eigen {
  template<class> struct NumTraits;
  template<> struct NumTraits<CORE::BigInt>
  {
    typedef CORE::BigInt Real;
    typedef CORE::BigRat NonInteger;
    typedef CORE::BigInt Nested;
    typedef CORE::BigInt Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = 1,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 30,
      MulCost = 50
    };
  };
}

#include <CGAL/enable_warnings.h>

#endif // CGAL_CORE_BIGINT_H
