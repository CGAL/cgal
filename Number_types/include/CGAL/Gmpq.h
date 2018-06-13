// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>

#ifndef CGAL_GMPQ_H
#define CGAL_GMPQ_H

#include <CGAL/number_type_basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmp_coercion_traits.h>

namespace CGAL {

// AST for Gmpq-class
template <> class Algebraic_structure_traits< Gmpq >
  : public Algebraic_structure_traits_base< Gmpq, Field_tag >  {
  public:
    typedef Tag_true            Is_exact;
    typedef Tag_false            Is_numerical_sensitive;

    class Is_square
      : public CGAL::binary_function< Type, Type&,
                                bool > {
      public:
        bool operator()( const Type& x_, Type& y ) const {
          Gmpq x( x_ );
          mpq_canonicalize( x.mpq() );
          Algebraic_structure_traits< Gmpz >::Sqrt sqrt;
          y = Gmpq( sqrt( x.numerator() ), sqrt( x.denominator() ) );
          return y*y == x;
        }
        bool operator()( const Type& x) const {
            Type y;
            return operator()(x,y);
        }

    };

    class Simplify
      : public CGAL::unary_function< Type&, void > {
      public:
        void operator()( Type& x) const {
          mpq_canonicalize( x.mpq() );
        }
    };

};

// RET for Gmpq-class

template <> class Real_embeddable_traits< Gmpq >
  : public INTERN_RET::Real_embeddable_traits_base< Gmpq , CGAL::Tag_true > {
  public:
  
    class Sgn
      : public CGAL::unary_function< Type, ::CGAL::Sign > {
      public:
        ::CGAL::Sign operator()( const Type& x ) const {
          return x.sign();
        }
    };

    class To_double
      : public CGAL::unary_function< Type, double > {
      public:
        double operator()( const Type& x ) const {
          return x.to_double();
        }
    };

    class To_interval
      : public CGAL::unary_function< Type, std::pair< double, double > > {
      public:
        std::pair<double, double> operator()( const Type& x ) const {
#if MPFR_VERSION_MAJOR >= 3
	  mpfr_exp_t emin = mpfr_get_emin();
	  mpfr_set_emin(-1073);
	  MPFR_DECL_INIT (y, 53); /* Assume IEEE-754 */
	  int r = mpfr_set_q (y, x.mpq(), MPFR_RNDA);
	  r = mpfr_subnormalize (y, r, MPFR_RNDA); /* Round subnormals */
	  double i = mpfr_get_d (y, MPFR_RNDA); /* EXACT but can overflow */
	  mpfr_set_emin(emin); /* Restore old value, users may care */
	  // With mpfr_set_emax(1024) we could drop the is_finite test
	  if (r == 0 && is_finite (i))
	    return std::pair<double, double>(i, i);
	  else
	    {
	      double s = nextafter (i, 0);
	      if (i < 0)
		return std::pair<double, double>(i, s);
	      else
		return std::pair<double, double>(s, i);
	    }
#else
          mpfr_t y;
          mpfr_init2 (y, 53); /* Assume IEEE-754 */
          mpfr_set_q (y, x.mpq(), GMP_RNDD);
          double i = mpfr_get_d (y, GMP_RNDD); /* EXACT but can overflow */
          mpfr_set_q (y, x.mpq(), GMP_RNDU);
          double s = mpfr_get_d (y, GMP_RNDU); /* EXACT but can overflow */
          mpfr_clear (y);
          return std::pair<double, double>(i, s);
#endif
        }
    };
};

/*! \ingroup NiX_Fraction_traits_spec
 *  \brief Specialization of Fraction_traits for Gmpq
 */
template <>
class Fraction_traits< Gmpq > {
public:
    typedef Gmpq Type;
    typedef ::CGAL::Tag_true Is_fraction;
    typedef Gmpz Numerator_type;
    typedef Gmpz Denominator_type;
    typedef Algebraic_structure_traits< Gmpz >::Gcd Common_factor;
    class Decompose {
    public:
        typedef Gmpq first_argument_type;
        typedef Gmpz& second_argument_type;
        typedef Gmpz& third_argument_type;
        void operator () (const Gmpq& rat, Gmpz& num,Gmpz& den) {
            num = rat.numerator();
            den = rat.denominator();
        }
    };
    class Compose {
    public:
        typedef Gmpz first_argument_type;
        typedef Gmpz second_argument_type;
        typedef Gmpq result_type;
        Gmpq operator () (const Gmpz& num,const Gmpz& den) {
            return Gmpq(num, den);
        }
    };
};

} //namespace CGAL

namespace Eigen {
  template<class> struct NumTraits;
  template<> struct NumTraits<CGAL::Gmpq>
  {
    typedef CGAL::Gmpq Real;
    typedef CGAL::Gmpq NonInteger;
    typedef CGAL::Gmpq Nested;
    typedef CGAL::Gmpq Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 150,
      MulCost = 100
    };
  };

  namespace internal {
    template<>
      struct significant_decimals_impl<CGAL::Gmpq>
      {
	static inline int run()
	{
	  return 0;
	}
      };
  }
}

//since types are included by Gmp_coercion_traits.h:
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpzf.h>
#include <CGAL/GMP_arithmetic_kernel.h>

#endif // CGAL_GMPQ_H
