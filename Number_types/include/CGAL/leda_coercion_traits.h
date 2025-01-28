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


/*! \file NiX/LEDA/Coercion_traits.h
 *  \brief provides specializations of Coercion_traits for the LEDA number types.
 */

#ifndef CGAL_LEDA_COERCION_TRAITS_H
#define CGAL_LEDA_COERCION_TRAITS_H

#include <CGAL/number_type_basic.h>

#ifdef CGAL_USE_LEDA

#include <CGAL/LEDA_basic.h>
#include <LEDA/numbers/integer.h>
#include <LEDA/numbers/bigfloat.h>
#include <LEDA/numbers/rational.h>
#include <LEDA/numbers/real.h>

namespace CGAL {


//LEDA internal coercions:

    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::integer,::leda::bigfloat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::integer,::leda::rational)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::integer,::leda::real)

// CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::bigfloat,::leda::rational); see below
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::bigfloat,::leda::real)

    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(::leda::rational,::leda::real)

// The following definitions reflect the interaction of the LEDA number types
// with the built in types,
// leda integer:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short    ,::leda::integer)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int      ,::leda::integer)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long     ,::leda::integer)

// leda rational:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,::leda::rational)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,::leda::rational)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,::leda::rational)

    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,::leda::rational)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,::leda::rational)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::leda::rational)

// leda bigfloat:      :
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,::leda::bigfloat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,::leda::bigfloat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,::leda::bigfloat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,::leda::bigfloat)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,::leda::bigfloat)

// leda real:
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(short      ,::leda::real)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(int        ,::leda::real)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(float      ,::leda::real)
    CGAL_DEFINE_COERCION_TRAITS_FROM_TO(double     ,::leda::real)


//not provided by LEDA
//Note that this is not symmetric to CORE
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long,::leda::integer);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long,::leda::rational);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::leda::bigfloat);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::leda::bigfloat);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long       ,::leda::real);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long long  ,::leda::real);
//CGAL_DEFINE_COERCION_TRAITS_FROM_TO(long double,::leda::real);


template <>
struct Coercion_traits< ::leda::bigfloat ,::leda::rational  >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_false Are_implicit_interoperable;
    typedef ::leda::rational Type;
    struct Cast{
        typedef Type result_type;
        Type operator()(const ::leda::rational& x)  const { return x;}
        Type operator()(const ::leda::bigfloat& x) const {
            return x.to_rational();
        }
    };
};
template <> struct Coercion_traits< ::leda::rational, ::leda::bigfloat >
    :public Coercion_traits< ::leda::bigfloat , ::leda::rational >{};


} //namespace CGAL
#endif // CGAL_USE_LEDA
#endif //CGAL_LEDA_COERCION_TRAITS_H
//EOF
