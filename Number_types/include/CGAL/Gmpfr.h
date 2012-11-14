// Copyright (c) 2007-2009 Inria Lorraine (France). All rights reserved.
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
        public std::unary_function<Type,Boolean>{
                Boolean operator()(const Type &x)const{
                        return x.is_zero();
                }
        };

        struct Is_one:
        public std::unary_function<Type,Boolean>{
                Boolean operator()(const Type &x)const{
                        return x.is_one();
                }
        };

        struct Square:
        public std::unary_function<Type,Type>{
                Type operator()(const Type &x)const{
                        return x.square();
                };
        };

        struct Is_square:
        public std::binary_function<Type,Type&,Boolean>{
                Boolean operator()(const Type &x,Type &y)const{
                        return x.is_square(y);
                };
                Boolean operator()(const Type &x)const{
                        return x.is_square();
                };
        };

        struct Sqrt:
        public std::unary_function<Type,Type>{
                Type operator()(const Type &x)const{
                        return x.sqrt();
                };
        };

        struct Kth_Root:
        public std::binary_function<int,Type,Type>{
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
        public std::unary_function<Type,Boolean>{
                inline Boolean operator()(const Type &x)const{
                        return(x.is_number());
                };
        };

        struct Abs:
        public std::unary_function<Type,Type>{
                inline Type operator()(const Type &x)const{
                        return x.abs();
                };
        };

        struct Sgn:
        public std::unary_function<Type,Sign>{
                inline Sign operator()(const Type &x)const{
                        return x.sign();
                };
        };

        struct Is_positive:
        public std::unary_function<Type,Boolean>{
                inline Boolean operator()(const Type &x)const{
                        return(x.sign()==POSITIVE);
                };
        };

        struct Is_negative:
        public std::unary_function<Type,Boolean>{
                inline Boolean operator()(const Type &x)const{
                        return(x.sign()==NEGATIVE);
                };
        };

        struct Compare:
        public std::binary_function<Type,Type,Comparison_result>{
                inline Comparison_result operator()
                        (const Type &x,const Type &y)const{
                                return x.compare(y);
                        }; 
          CGAL_IMPLICIT_INTEROPERABLE_BINARY_OPERATOR_WITH_RT(Type,Comparison_result)
        };

        struct To_double:
        public std::unary_function<Type,double>{
                inline double operator()(const Type &x)const{
                        return x.to_double();
                };
        };

        struct To_interval:
        public std::unary_function<Type,std::pair<double,double> >{
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

    static inline Real epsilon() { return 0; }

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
