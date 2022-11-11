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


#ifndef CGAL_CORE_COERCION_TRAITS_H
#define CGAL_CORE_COERCION_TRAITS_H 1

#include <CGAL/number_type_basic.h>

#ifdef CGAL_USE_CORE

#ifndef CORE_LEVEL
#define CORE_LEVEL 4
#endif

#include <CGAL/CORE/CORE.h>

//#include <NiX/Coercion_traits.h>

namespace CGAL {

//CORE internal coercions:


// The following definitions reflect the interaction of the CORE number types
// with the built in types,
// CORE BigInt:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short    ,::CORE::BigInt)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,::CORE::BigInt)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long     ,::CORE::BigInt)


// CORE BigRat:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short         ,::CORE::BigRat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int           ,::CORE::BigRat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long          ,::CORE::BigRat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float         ,::CORE::BigRat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double        ,::CORE::BigRat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::CORE::BigInt,::CORE::BigRat)

// CORE Expr:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short         ,::CORE::Expr)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int           ,::CORE::Expr)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long          ,::CORE::Expr)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float         ,::CORE::Expr)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double        ,::CORE::Expr)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::CORE::BigInt,::CORE::Expr)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::CORE::BigRat,::CORE::Expr)



// NOTE that CORE::BigFloat counts as an interval !
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short         ,::CORE::BigFloat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int           ,::CORE::BigFloat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long          ,::CORE::BigFloat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float         ,::CORE::BigFloat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double        ,::CORE::BigFloat)

template <>
struct Coercion_traits<CORE::BigFloat , ::CORE::BigInt>{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    typedef CORE::BigFloat Type;

    struct Cast{
        typedef Type result_type;
        Type operator()(const CORE::BigFloat& x)  const { return x;}
        Type operator()(const ::CORE::BigInt x) const {
            CORE::BigFloat result;
            result.approx(x,CORE::get_static_defRelPrec().toLong(),LONG_MAX);
            // Do not use MakeFloorExact as it changes the Bigfloat
            CGAL_postcondition( ::CORE::BigRat(::CORE::BigFloat(result.m()-result.err(),0,result.exp())) <= x );
            CGAL_postcondition( ::CORE::BigRat(::CORE::BigFloat(result.m()+result.err(),0,result.exp())) >= x );
            return result;
        }
    };
};

template <>
struct Coercion_traits<CORE::BigFloat , ::CORE::BigRat>{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    typedef CORE::BigFloat Type;

    struct Cast{
        typedef Type result_type;
        Type operator()(const CORE::BigFloat& x)  const { return x;}
        Type operator()(const ::CORE::BigRat x) const {

          CORE::BigFloat result(x,CORE::get_static_defRelPrec().toLong(),LONG_MAX);
            // Do not use MakeFloorExact as it changes the Bigfloat
            CGAL_postcondition( ::CORE::BigRat(::CORE::BigFloat(result.m()-result.err(),0,result.exp())) <= x );
            CGAL_postcondition( ::CORE::BigRat(::CORE::BigFloat(result.m()+result.err(),0,result.exp())) >= x );
            return result;
        }
    };
};

template <>
struct Coercion_traits<CORE::BigFloat , ::CORE::Expr>{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    typedef CORE::BigFloat Type;

    struct Cast{
        typedef Type result_type;
        Type operator()(const CORE::BigFloat& x)  const { return x;}
        Type operator()(const ::CORE::Expr x) const {
            CORE::BigFloat result(x, CORE::get_static_defRelPrec().toLong(),LONG_MAX);
            // Do not use MakeFloorExact as it changes the Bigfloat
            CGAL_postcondition( ::CORE::BigRat(::CORE::BigFloat(result.m()-result.err(),0,result.exp())) <= x );
            CGAL_postcondition( ::CORE::BigRat(::CORE::BigFloat(result.m()+result.err(),0,result.exp())) >= x );
            return result;
        }
    };
};

template <> struct Coercion_traits< ::CORE::BigInt, CORE::BigFloat >
    :public Coercion_traits< CORE::BigFloat , ::CORE::BigInt >{};

template <> struct Coercion_traits< ::CORE::BigRat, CORE::BigFloat >
    :public Coercion_traits< CORE::BigFloat , ::CORE::BigRat >{};

template <> struct Coercion_traits< ::CORE::Expr, CORE::BigFloat >
    :public Coercion_traits< CORE::BigFloat , ::CORE::Expr>{};






// not provieded by CORE
// Note that this is not symmetric to LEDA
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::CORE::BigInt)
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::CORE::BigRat)
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::CORE::BigRat)
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::CORE::BigFloat)
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::CORE::BigFloat)
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::CORE::Expr)
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::CORE::Expr)

} //namespace CGAL

#endif // CGAL_USE_CORE
#endif //CGAL_CORE_COERCION_TRAITS_H 1
//EOF
