//
//  Copyright (c) 2000-2002
//  Joerg Walter, Mathias Koch
//
//  Permission to use, copy, modify, distribute and sell this software
//  and its documentation for any purpose is hereby granted without fee,
//  provided that the above copyright notice appear in all copies and
//  that both that copyright notice and this permission notice appear
//  in supporting documentation.  The authors make no representations
//  about the suitability of this software for any purpose.
//  It is provided "as is" without express or implied warranty.
//
//  The authors gratefully acknowledge the support of
//  GeNeSys mbH & Co. KG in producing this work.
//

#ifndef BOOST_UBLAS_CONFIG_H
#define BOOST_UBLAS_CONFIG_H

#include <cassert>
#include <cstddef>

#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/limits.hpp>
#include <boost/utility.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_reference.hpp>

// Compiler specific problems: default configuration
#if defined (BOOST_STRICT_CONFIG) || ! (\
        defined (BOOST_MSVC) || \
        defined (__GNUC__) || \
        defined (__BORLANDC__) || \
        defined (_ICL) || \
        defined (_ICC) || \
        defined (__COMO__) || \
        defined (__MWERKS__))

#define BOOST_UBLAS_TYPENAME typename
#define BOOST_UBLAS_USING using
// This could be eliminated.
#define BOOST_UBLAS_EXPLICIT explicit

#define BOOST_UBLAS_USE_LONG_DOUBLE

#define BOOST_UBLAS_USE_STREAM

#endif



// Microsoft Visual C++
#if defined (BOOST_MSVC) && ! defined (BOOST_STRICT_CONFIG)

// Version 6.0 & 7.0
#if BOOST_MSVC <= 1300

// Disable some MSVC specific warnings.
#pragma warning (disable: 4355)
#pragma warning (disable: 4503)
#pragma warning (disable: 4786)

// MSVC doesn't always accept the 'typename' keyword.
#define BOOST_UBLAS_TYPENAME
// MSVC doesn't accept the 'using' keyword (at least for importing base members).
#define BOOST_UBLAS_USING
#define BOOST_UBLAS_EXPLICIT explicit

// With MSVC we could perform IO via basic_stream
// #define BOOST_UBLAS_USE_BASIC_STREAM
// IO via streams
#define BOOST_UBLAS_USE_STREAM

// MSVC extensions seem to disable abs () overloads in <cmath>.
#ifdef _MSC_EXTENSIONS
#define BOOST_UBLAS_NO_CMATH
#endif

#define BOOST_UBLAS_NO_ELEMENT_PROXIES
// This seems to work now thanks to the great work of the MPL team!
// #define BOOST_UBLAS_NO_SMART_PROXIES

// Using MSVC the following is missing:
// namespace std {
//
//     typedef unsigned size_t;
//     typedef signed ptrdiff_t;
//
//     template<class T>
//     const T &max (const T &t1, const T &t2) {
//         return t1 > t2 ? t1 : t2;
//     }
//     template<class T>
//     const T &min (const T &t1, const T &t2) {
//         return t1 < t2 ? t1 : t2;
//     }
//
// }
// But boost provides it for us.

// This seems to be a problem in boost.config, but won't be fixed.
#ifdef __SGI_STL_PORT
#define BOOST_MSVC_STD_ITERATOR
#endif

// Version 7.1
#else

#define BOOST_UBLAS_TYPENAME typename
#define BOOST_UBLAS_USING using
// This could be eliminated.
#define BOOST_UBLAS_EXPLICIT explicit

#define BOOST_UBLAS_USE_LONG_DOUBLE

#define BOOST_UBLAS_USE_STREAM

// One of these workarounds is needed for MSVC 7.1 AFAIK
// (thanks to John Maddock and Martin Lauer).
// The second workaround looks like BOOST_UBLAS_QUALIFIED_TYPENAME.
#define BOOST_UBLAS_NO_NESTED_CLASS_RELATION
// #define BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION

#endif

#endif



// GNU Compiler Collection
#if defined (__GNUC__) && ! defined (BOOST_STRICT_CONFIG)

#define BOOST_UBLAS_TYPENAME typename
#define BOOST_UBLAS_USING using
#define BOOST_UBLAS_EXPLICIT explicit

// GCC 2.95.3 allows to use iterator_base_traits.
#define BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS
// GCC 2.95.3 needs BOOST_UBLAS_REVERSE_ITERATOR_OVERLOADS (this seems to be arguable).
#define BOOST_UBLAS_REVERSE_ITERATOR_OVERLOADS

#define BOOST_UBLAS_USE_LONG_DOUBLE

#define BOOST_UBLAS_USE_STREAM

#if __GNUC__ <= 2 && __GNUC_MINOR__ <= 95
#define BOOST_UBLAS_NO_MEMBER_FRIENDS
#define BOOST_UBLAS_NO_PROXY_SHORTCUTS
#endif

#endif



// Borland Compiler
#if defined (__BORLANDC__) && ! defined (BOOST_STRICT_CONFIG)

#define BOOST_UBLAS_TYPENAME typename
#define BOOST_UBLAS_USING using
#define BOOST_UBLAS_EXPLICIT explicit

#define BOOST_UBLAS_USE_LONG_DOUBLE

#define BOOST_UBLAS_USE_STREAM

#define BOOST_UBLAS_NO_ELEMENT_PROXIES
#define BOOST_UBLAS_NO_SMART_PROXIES

#define BOOST_UBLAS_NO_PROXY_SHORTCUTS
#define BOOST_UBLAS_NO_DERIVED_HELPERS

// BCC's <complex> broken.
// Thanks to John Maddock for providing a workaround.
#if defined(_STLPORT_VERSION) && defined(_STLP_USE_OWN_NAMESPACE) && !defined(std)
#include <complex>
namespace std {
    using stlport::abs;
    using stlport::sqrt;
}
#endif

#endif



// Thanks to Roberto Andres Ruiz Vial for porting to Intel.
#if defined (__ICL) && ! defined (BOOST_STRICT_CONFIG)

#define BOOST_UBLAS_TYPENAME typename
#define BOOST_UBLAS_USING using
#define BOOST_UBLAS_EXPLICIT explicit

#define BOOST_UBLAS_USE_LONG_DOUBLE

#define BOOST_UBLAS_USE_STREAM

#endif



// Intel Compiler under Linux
#if defined (__ICC) && ! defined (BOOST_STRICT_CONFIG)

#define BOOST_UBLAS_TYPENAME typename
#define BOOST_UBLAS_USING using
#define BOOST_UBLAS_EXPLICIT explicit

// ICC sometimes needs qualified type names.
#define BOOST_UBLAS_QUALIFIED_TYPENAME

#define BOOST_UBLAS_USE_LONG_DOUBLE

#define BOOST_UBLAS_USE_STREAM

#define BOOST_UBLAS_USE_SIMD

namespace boost { namespace numeric { namespace ublas {

    template<class C, class IC>
    class indexed_iterator;

    template<class V>
    class index_pair;
    template<class M>
    class index_triple;

}}}

namespace std {

    // Needed for Intel on Itanium?
    template<class C, class IC>
    inline
    void iter_swap (boost::numeric::ublas::indexed_iterator<C, IC> it1,
                    boost::numeric::ublas::indexed_iterator<C, IC> it2) {
        swap (*it1, *it2);
    }

    template<class V>
    inline
    void swap (boost::numeric::ublas::index_pair<V> i1, boost::numeric::ublas::index_pair<V> i2) {
        i1.swap (i2);
    }
    template<class M>
    inline
    void swap (boost::numeric::ublas::index_triple<M> i1, boost::numeric::ublas::index_triple<M> i2) {
        i1.swap (i2);
    }

}

#endif



// Thanks to Kresimir Fresl for porting to Comeau.
#if defined (__COMO__) && ! defined (BOOST_STRICT_CONFIG)

#define BOOST_UBLAS_TYPENAME typename
#define BOOST_UBLAS_USING using
#define BOOST_UBLAS_EXPLICIT explicit

// Comeau allows to use iterator_base_traits.
#define BOOST_UBLAS_USE_ITERATOR_BASE_TRAITS

#define BOOST_UBLAS_USE_LONG_DOUBLE

#define BOOST_UBLAS_USE_STREAM

#endif



// Metrowerks Codewarrior
#if defined (__MWERKS__) && ! defined (BOOST_STRICT_CONFIG)

#define BOOST_UBLAS_TYPENAME typename
#define BOOST_UBLAS_USING using
#define BOOST_UBLAS_EXPLICIT explicit

#define BOOST_UBLAS_USE_LONG_DOUBLE

#define BOOST_UBLAS_USE_STREAM

#define BOOST_UBLAS_NO_MEMBER_FRIENDS

#endif



// Enable assignment of non conformant proxies
// Thanks to Michael Stevens for spotting this.
#define BOOST_UBLAS_NON_CONFORMANT_PROXIES

// Enable different sparse element proxies
// These fix a [1] = a [0] = 1, but probably won't work on broken compilers.
// Thanks to Marc Duflot for spotting this.
#ifndef BOOST_UBLAS_NO_ELEMENT_PROXIES
// #define BOOST_UBLAS_STRICT_STORAGE_SPARSE
#define BOOST_UBLAS_STRICT_VECTOR_SPARSE
#define BOOST_UBLAS_STRICT_MATRIX_SPARSE
// #define BOOST_UBLAS_STRICT_HERMITIAN
#endif

// Enable compile time typedefs for proxies
#ifndef BOOST_UBLAS_NO_SMART_PROXIES
#define BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
// #define BOOST_UBLAS_REFERENCE_CONST_MEMBER
#define BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
// #define BOOST_UBLAS_PROXY_CONST_MEMBER
#define BOOST_UBLAS_CT_PROXY_CLOSURE_TYPEDEFS
#endif



// Enable performance options in release mode
#ifdef NDEBUG

#ifdef BOOST_MSVC
// MSVC has special inlining options
#pragma inline_recursion (on)
#pragma inline_depth (255)
#pragma auto_inline (on)
// #define BOOST_UBLAS_INLINE __forceinline
#define BOOST_UBLAS_INLINE __inline
#else
#define BOOST_UBLAS_INLINE inline
#endif

// Do not check sizes!
#define BOOST_UBLAS_USE_FAST_SAME

// Use expression templates.
#define BOOST_UBLAS_USE_ET

// Disable performance options in debug mode
#else

#ifdef BOOST_MSVC
// MSVC has special inlining options
// #pragma inline_recursion (off)
// #pragma inline_depth ()
// #pragma auto_inline (off)
#endif
#define BOOST_UBLAS_INLINE

#ifdef BOOST_MSVC
// Use expression templates (otherwise we get many ICE's)
#define BOOST_UBLAS_USE_ET
#endif

// Bounds check
#define BOOST_UBLAS_BOUNDS_CHECK
// Type check for non dense matrices
#define BOOST_UBLAS_TYPE_CHECK

#endif

#ifdef BOOST_UBLAS_TYPE_CHECK
static bool disable_type_check = false;
#endif
#ifndef BOOST_UBLAS_TYPE_CHECK_EPSILON
#define BOOST_UBLAS_TYPE_CHECK_EPSILON (type_traits<real_type>::sqrt (std::numeric_limits<real_type>::epsilon ()))
#endif
#ifndef BOOST_UBLAS_TYPE_CHECK_MIN
#define BOOST_UBLAS_TYPE_CHECK_MIN (type_traits<real_type>::sqrt (std::numeric_limits<real_type>::min ()))
#endif



// Use invariant hoisting.
// #define BOOST_UBLAS_USE_INVARIANT_HOISTING

// Use Duff's device
// #define BOOST_UBLAS_USE_DUFF_DEVICE

// Choose evaluation method for dense vectors and matrices
#define BOOST_UBLAS_USE_INDEXING
// #define BOOST_UBLAS_USE_ITERATING
// #define BOOST_UBLAS_ITERATOR_THRESHOLD 0
// #define BOOST_UBLAS_ITERATOR_THRESHOLD (std::numeric_limits<std::ptrdiff_t>::max ())

// Use indexed iterators.
// #define BOOST_UBLAS_USE_INDEXED_ITERATOR

// 16 byte aligned arrays (for ICC)
// #define BOOST_UBLAS_ALIGN_16 __declspec (align (16))
#define BOOST_UBLAS_ALIGN_16

#include <boost/numeric/ublas/fwd.hpp>

#endif




























