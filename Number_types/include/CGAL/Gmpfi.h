// Copyright (c) 2007-2009 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>
//         Michael Hemmer <Michael.Hemmer@sophia.inria.fr>

#ifndef CGAL_GMPFI_H
#define CGAL_GMPFI_H

#include <CGAL/config.h>
#ifdef CGAL_USE_MPFI

#include <CGAL/GMP/Gmpfi_type.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/mpfi_coercion_traits.h>
#include <CGAL/Interval_traits.h>
#include <CGAL/Bigfloat_interval_traits.h>

namespace CGAL{

template <>
class Algebraic_structure_traits<Gmpfi>:
public Algebraic_structure_traits_base<Gmpfi,Field_with_kth_root_tag>{
public:

        typedef Tag_false       Is_exact;
        typedef Tag_true        Is_numerical_sensitive;
        typedef Uncertain<bool> Boolean;

        struct Is_zero:
        public CGAL::cpp98::unary_function<Type,Boolean>{
                Boolean operator()(const Type &x)const{
                        return x.is_zero();
                }
        };

        struct Is_one:
        public CGAL::cpp98::unary_function<Type,Boolean>{
                Boolean operator()(const Type &x)const{
                        return x.is_one();
                }
        };

        struct Square:
        public CGAL::cpp98::unary_function<Type,Type>{
                Type operator()(const Type &x)const{
                        return x.square();
                };
        };

        struct Is_square:
        public CGAL::cpp98::binary_function<Type,Type&,Boolean>{
                Boolean operator()(const Type &x)const{
                        return x.is_square();
                };
                Boolean operator()(const Type &x,Type &y)const{
                        return x.is_square(y);
                };
        };

        struct Sqrt:
        public CGAL::cpp98::unary_function<Type,Type>{
                Type operator()(const Type &x)const{
                        return x.sqrt();
                };
        };

        struct Kth_Root:
        public CGAL::cpp98::binary_function<int,Type,Type>{
                Type operator()(int k,const Type &x)const{
                        return (k==3?x.cbrt():x.kthroot(k));
                };
        };

        struct Divides:
        public CGAL::cpp98::binary_function<Type,Type,Boolean>{
                Boolean operator()(const Type &d,const Type &n)const{
                        // Avoid compiler warning
                        (void)n;
                        return !(d.is_zero());
                };
                Boolean operator()(const Type &d,const Type &n,Type &c)const{
                        return d.divides(n,c);
                };
        };
};

template <>
class Real_embeddable_traits<Gmpfi>:
public INTERN_RET::Real_embeddable_traits_base<Gmpfi,CGAL::Tag_true>{

        typedef Algebraic_structure_traits<Type>        AST;

        public:

        typedef Tag_true                                Is_real_embeddable;
        typedef Uncertain<bool>                         Boolean;
        typedef Uncertain<CGAL::Comparison_result>      Comparison_result;
        typedef Uncertain<CGAL::Sign>                   Sign;

        typedef AST::Is_zero    Is_zero;

        struct Is_finite:
        public CGAL::cpp98::unary_function<Type,Boolean>{
                inline Boolean operator()(const Type &x)const{
                        return(x.is_number());
                };
        };

        struct Abs:
        public CGAL::cpp98::unary_function<Type,Type>{
                inline Type operator()(const Type &x)const{
                        return x.abs();
                };
        };

        struct Sgn:
        public CGAL::cpp98::unary_function<Type,Sign>{
                inline Sign operator()(const Type &x)const{
                        return x.sign();
                };
        };

        struct Is_positive:
        public CGAL::cpp98::unary_function<Type,Boolean>{
                inline Boolean operator()(const Type &x)const{
                        return x.is_positive();
                };
        };

        struct Is_negative:
        public CGAL::cpp98::unary_function<Type,Boolean>{
                inline Boolean operator()(const Type &x)const{
                        return x.is_negative();
                };
        };

        struct Compare:
        public CGAL::cpp98::binary_function<Type,Type,Comparison_result>{
                inline Comparison_result operator()
                        (const Type &x,const Type &y)const{
                                return x.compare(y);
                        };
          CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT(Type,Comparison_result)
        };

        struct To_double:
        public CGAL::cpp98::unary_function<Type,double>{
                inline double operator()(const Type &x)const{
                        return x.to_double();
                };
        };

        struct To_interval:
        public CGAL::cpp98::unary_function<Type,std::pair<double,double> >{
                inline std::pair<double,double> operator()(const Type &x)const{
                                return x.to_interval();
                        };
                };
};



template<>
class Interval_traits<Gmpfi>
  : public internal::Interval_traits_base<Gmpfi>{
public:
  typedef Interval_traits<Gmpfi> Self;
  typedef Gmpfi Interval;
  typedef CGAL::Gmpfr Bound;
  typedef CGAL::Tag_false With_empty_interval;
  typedef CGAL::Tag_true  Is_interval;

  struct Construct :public CGAL::cpp98::binary_function<Bound,Bound,Interval>{
    Interval operator()( const Bound& l,const Bound& r) const {
      CGAL_precondition( l < r );
      return Interval(std::make_pair(l,r));
    }
  };

  struct Lower :public CGAL::cpp98::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.inf();
    }
  };

  struct Upper :public CGAL::cpp98::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.sup();
    }
  };

  struct Width :public CGAL::cpp98::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return Gmpfr::sub(a.sup(),a.inf(),std::round_toward_infinity);
    }
  };

  struct Median :public CGAL::cpp98::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return (a.inf()+a.sup())/2;
    }
  };

  struct Norm :public CGAL::cpp98::unary_function<Interval,Bound>{
    Bound operator()( const Interval& a ) const {
      return a.abs().sup();
    }
  };

  struct Singleton :public CGAL::cpp98::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return a.inf() == a.sup();
    }
  };

  struct Zero_in :public CGAL::cpp98::unary_function<Interval,bool>{
    bool operator()( const Interval& a ) const {
      return a.inf() <= 0  &&  0 <= a.sup();
    }
  };

  struct In :public CGAL::cpp98::binary_function<Bound,Interval,bool>{
    bool operator()( Bound x, const Interval& a ) const {
      return a.inf() <= x && x <= a.sup();
    }
  };

  struct Equal :public CGAL::cpp98::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return a.is_same(b);
    }
  };

  struct Overlap :public CGAL::cpp98::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return a.do_overlap(b);
    }
  };

  struct Subset :public CGAL::cpp98::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return b.inf() <= a.inf() && a.sup() <= b.sup() ;
    }
  };

  struct Proper_subset :public CGAL::cpp98::binary_function<Interval,Interval,bool>{
    bool operator()( const Interval& a, const Interval& b ) const {
      return Subset()(a,b) && ! Equal()(a,b);
    }
  };

  struct Hull :public CGAL::cpp98::binary_function<Interval,Interval,Interval>{
    Interval operator()( const Interval& a, const Interval& b ) const {
      BOOST_USING_STD_MAX();
      BOOST_USING_STD_MIN();
      return Interval(
          std::make_pair(
              min BOOST_PREVENT_MACRO_SUBSTITUTION (a.inf(),b.inf()),
              max BOOST_PREVENT_MACRO_SUBSTITUTION (a.sup(),b.sup())));
    }
  };


//  struct Empty is Null_functor

  struct Intersection :public CGAL::cpp98::binary_function<Interval,Interval,Interval>{
    Interval operator()( const Interval& a, const Interval& b ) const {
      BOOST_USING_STD_MAX();
      BOOST_USING_STD_MIN();
      Bound l(max BOOST_PREVENT_MACRO_SUBSTITUTION (Lower()(a),Lower()(b)));
      Bound u(min BOOST_PREVENT_MACRO_SUBSTITUTION (Upper()(a),Upper()(b)));
      if(u < l ) throw Exception_intersection_is_empty();
      return Construct()(l,u);
    }
  };
};

template<>
class Bigfloat_interval_traits<Gmpfi>
  :public Interval_traits<Gmpfi>
{
  typedef Gmpfi NT;
  typedef CGAL::Gmpfr BF;
public:
  typedef Bigfloat_interval_traits<Gmpfi> Self;
  typedef CGAL::Tag_true                  Is_bigfloat_interval;

  struct Relative_precision: public CGAL::cpp98::unary_function<NT,long>{

    long operator()(const NT& x) const {
      CGAL_precondition(!Singleton()(x));
      CGAL_precondition(!CGAL::zero_in(x));

      // w = |x| * 2^-p (return p)
      BF w(CGAL::width(x));
      mpfr_div(w.fr(), w.fr(), CGAL::lower(CGAL::abs(x)).fr(), GMP_RNDU);
      mpfr_log2(w.fr(), w.fr(), GMP_RNDU);
      return -mpfr_get_si(w.fr(), GMP_RNDU);
    }
  };

  struct Set_precision {
    // type for the \c AdaptableUnaryFunction concept.
    typedef long  argument_type;
    // type for the \c AdaptableUnaryFunction concept.
    typedef long  result_type;

    long operator()( long prec ) const {
      return Gmpfi::set_default_precision(prec);
    }
  };

  struct Get_precision {
    // type for the \c AdaptableGenerator concept.
    typedef long  result_type;
    long operator()() const {
      return Gmpfi::get_default_precision();
    }
  };
};

} // namespace CGAL

namespace Eigen {
  template<class> struct NumTraits;
  template<> struct NumTraits<CGAL::Gmpfi>
  {
    typedef CGAL::Gmpfi Real;
    typedef CGAL::Gmpfi NonInteger;
    typedef CGAL::Gmpfi Nested;
    typedef CGAL::Gmpfi Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 12,
      AddCost = 100,
      MulCost = 100
    };
  };
}

#include <CGAL/GMP/Gmpfi_type.h>
#include <CGAL/GMP_arithmetic_kernel.h>

#endif

#endif  // CGAL_GMPFI_H

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab
