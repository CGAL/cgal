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
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>



/*! \file NiX/Gmp/Coercion_traits.h
 *  \brief Provides specializations of Coercion_traits for the Gmp types.
 */

#ifndef CGAL_GMPXX_COERCION_TRAITS_H
#define CGAL_GMPXX_COERCION_TRAITS_H 1

#include <CGAL/number_type_basic.h>
#include <CGAL/Coercion_traits.h>

#include <cstring> // needed by GMP 4.1.4 since <gmpxx.h> misses it.
#include <gmpxx.h>
#include <mpfr.h>

namespace CGAL {

//mpz_class internal coercions:
//self for mpz_class / mpq_class
template <class T , class U>
struct Coercion_traits<
  ::__gmp_expr< T , U>,::__gmp_expr< T , U>  >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef ::__gmp_expr<T , T> Type;
    struct Cast{
        typedef Type result_type;
        template <class U3>
        Type operator()(const ::__gmp_expr< T , U3>& x) const {
            return x;
        }
    };
};

template <class T, class U1, class U2>
struct Coercion_traits<
  ::__gmp_expr< T , U1>,::__gmp_expr< T , U2>  >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef ::__gmp_expr< T , T > Type;
    struct Cast{
        typedef Type result_type;
        template <class U3>
        Type operator()(const ::__gmp_expr< T , U3>& x) const {
            return x;
        }
    };
};


template <class T1 , class T2, class U1, class U2>
struct Coercion_traits< ::__gmp_expr< T1 , U1>,::__gmp_expr< T2 , U2>  >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef mpq_class Type;
    struct Cast{
        typedef Type result_type;
        template <class T , class U>
        Type operator()(const ::__gmp_expr< T , U>& x) const {
            return Type(x);
        }
    };
};


// gmpzq_class implicit interoperable with int
template <class T, class U>
struct Coercion_traits<
  ::__gmp_expr< T , U >, int >{
    typedef Tag_true  Are_explicit_interoperable;
    typedef Tag_true  Are_implicit_interoperable;
    typedef ::__gmp_expr< T , T > Type;
    struct Cast{
        typedef Type result_type;
        template <class U3>
        Type operator()(const ::__gmp_expr< T , U3>& x) const {
            return x;
        }
        Type operator()(int x) const { return Type(x); }
    };
};
// gmpz_class implicit interoperable with int
template <class U, class T>
struct Coercion_traits< int , ::__gmp_expr< T , U> >
    :public Coercion_traits< ::__gmp_expr< T , U>, int >{};

} //namespace CGAL

#endif //CGAL_GMPXX_COERCION_TRAITS_H 1
//EOF
