// Copyright (c) 2007-2009 Inria Lorraine (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_GMPFR_H
#define CGAL_GMPFR_H

#include <CGAL/number_type_basic.h>
#include <CGAL/mpfr_coercion_traits.h>

namespace CGAL{

template <>
class Algebraic_structure_traits<Gmpfr>:
public Algebraic_structure_traits_base<Gmpfr,Field_with_kth_root_tag>{
public:

        typedef Tag_false       Is_exact;
        typedef Tag_true        Is_numerical_sensitive;
        typedef bool            Boolean;

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
                Boolean operator()(const Type &x,Type &y)const{
                        return x.is_square(y);
                };
                Boolean operator()(const Type &x)const{
                        return x.is_square();
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
};

template <>
class Real_embeddable_traits<Gmpfr>:
public INTERN_RET::Real_embeddable_traits_base<Gmpfr,CGAL::Tag_true>{

        typedef Algebraic_structure_traits<Type>        AST;

        public:

        typedef Tag_true                Is_real_embeddable;
        typedef bool                    Boolean;
        typedef CGAL::Comparison_result Comparison_result;
        typedef CGAL::Sign              Sign;

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
                        return(x.sign()==POSITIVE);
                };
        };

        struct Is_negative:
        public CGAL::cpp98::unary_function<Type,Boolean>{
                inline Boolean operator()(const Type &x)const{
                        return(x.sign()==NEGATIVE);
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
                inline std::pair<double,double>operator()(const Type &x)const{
                                return x.to_interval();
                };
        };
};

}

namespace Eigen {
  template<class> struct NumTraits;
  template<> struct NumTraits<CGAL::Gmpfr>
  {
    typedef CGAL::Gmpfr Real;
    typedef CGAL::Gmpfr NonInteger;
    typedef CGAL::Gmpfr Nested;
    typedef CGAL::Gmpfr Literal;

    static inline Real epsilon() { return 0; }
    static inline Real dummy_precision() { return 0; }

    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 6,
      AddCost = 50,
      MulCost = 50
    };
  };
}

#include <CGAL/GMP/Gmpfr_type.h>
#include <CGAL/GMP_arithmetic_kernel.h>

#endif  // CGAL_GMPFR_H

// vim: tabstop=8: softtabstop=8: smarttab: shiftwidth=8: expandtab
