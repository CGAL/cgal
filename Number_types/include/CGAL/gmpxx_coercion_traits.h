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



/*! \file NiX/Gmp/Coercion_traits.h
 *  \brief Provides specializations of Coercion_traits for the Gmp types.
 */

#ifndef CGAL_GMPXX_COERCION_TRAITS_H
#define CGAL_GMPXX_COERCION_TRAITS_H 1

#include <CGAL/Coercion_traits.h>

#include <cstring> // needed by GMP 4.1.4 since <gmpxx.h> misses it.
#include <gmpxx.h>
#include <mpfr.h>

namespace CGAL {

//self for mpz_class / mpq_class
template <class T, class U>
struct Coercion_traits< ::__gmp_expr< T , U>, ::__gmp_expr< T , U>  >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef ::__gmp_expr< T , T > Type;
    struct Cast{
        typedef Type result_type;
        template <class X>
        Type operator()(const X& x) const {
            return x;
        }
    };
};
template <class T, class U1, class U2>
struct Coercion_traits< ::__gmp_expr< T , U1>, ::__gmp_expr< T , U2>  >
: Coercion_traits< ::__gmp_expr<T,T>, ::__gmp_expr<T,T> > {};

//mixed mpz_class + mpq_class, ignore the possibility of mpf_class
template <class T1 , class T2, class U1, class U2>
struct Coercion_traits< ::__gmp_expr< T1 , U1>, ::__gmp_expr< T2 , U2>  >
: Coercion_traits< mpq_class, mpq_class > {};

// gmpzq_class implicit interoperable with int, short, long
template <class T, class U>
struct Coercion_traits< ::__gmp_expr< T , U >, int >
: public Coercion_traits< ::__gmp_expr< T , T >, ::__gmp_expr< T , T > > {};
template <class U, class T>
struct Coercion_traits< int , ::__gmp_expr< T , U> >
: public Coercion_traits< ::__gmp_expr< T , T >, ::__gmp_expr< T , T > > {};
template <class T, class U>
struct Coercion_traits< ::__gmp_expr< T , U >, long >
: public Coercion_traits< ::__gmp_expr< T , T >, ::__gmp_expr< T , T > > {};
template <class T, class U>
struct Coercion_traits< long , ::__gmp_expr< T , U > >
: public Coercion_traits< ::__gmp_expr< T , T >, ::__gmp_expr< T , T > > {};
template <class T, class U>
struct Coercion_traits< ::__gmp_expr< T , U >, short >
: public Coercion_traits< ::__gmp_expr< T , T >, ::__gmp_expr< T , T > > {};
template <class T, class U>
struct Coercion_traits< short , ::__gmp_expr< T , U > >
: public Coercion_traits< ::__gmp_expr< T , T >, ::__gmp_expr< T , T > > {};

// The traits are identical for float/double. The implicit conversion from double to mpz_class might disappear some day, but hopefully by that time CGAL will not support antediluvian versions of GMP anymore, which will make it easier to specialize the traits (or we will have a default version of Coercion_traits based on std::common_type and we can remove this file).
template <class T, class U>
struct Coercion_traits< ::__gmp_expr< T , U >, double >
: public Coercion_traits< ::__gmp_expr< T , T >, ::__gmp_expr< T , T > > {};
template <class T, class U>
struct Coercion_traits< double , ::__gmp_expr< T , U > >
: public Coercion_traits< ::__gmp_expr< T , T >, ::__gmp_expr< T , T > > {};
template <class T, class U>
struct Coercion_traits< ::__gmp_expr< T , U >, float >
: public Coercion_traits< ::__gmp_expr< T , T >, ::__gmp_expr< T , T > > {};
template <class T, class U>
struct Coercion_traits< float , ::__gmp_expr< T , U > >
: public Coercion_traits< ::__gmp_expr< T , T >, ::__gmp_expr< T , T > > {};

} //namespace CGAL

#endif //CGAL_GMPXX_COERCION_TRAITS_H 1
//EOF
