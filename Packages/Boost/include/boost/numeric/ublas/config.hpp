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

#include <memory>
#include <cassert>
#include <cstddef>

#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/limits.hpp>
#include <boost/noncopyable.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_reference.hpp>


// Microsoft Visual C++
#if defined (BOOST_MSVC) && ! defined (BOOST_STRICT_CONFIG)

// Version 6.0
#if BOOST_MSVC < 1300
// Cannot access private members from member class
#define BOOST_UBLAS_NESTED_CLASS_DR45
#endif

// Version 6.0 & 7.0
#if BOOST_MSVC <= 1300

// Disable some MSVC specific warnings.
#pragma warning (disable: 4355)
#pragma warning (disable: 4503)
#pragma warning (disable: 4786)

// static member data initialisers require = syntax
#define BOOST_UBLAS_STATIC_OLD_INIT
// Member friend syntax works but lookup sometime fails on VC7
#define BOOST_UBLAS_NO_MEMBER_FRIENDS
// Base traits templates syntax not supported
#define BOOST_UBLAS_NO_ITERATOR_BASE_TRAITS
// MSVC doesn't always accept the 'typename' keyword
#define BOOST_UBLAS_TYPENAME
// MSVC doesn't accept the 'using' keyword (at least for importing base members)
#define BOOST_UBLAS_USING
// MSVC doesn't support long double
#define BOOST_UBLAS_NO_LONG_DOUBLE

#ifdef NDEBUG
// MSVC has special inlining options
#pragma inline_recursion (on)
#pragma inline_depth (255)
#pragma auto_inline (on)
// #define BOOST_UBLAS_INLINE __forceinline
#define BOOST_UBLAS_INLINE __inline
#endif

// MSVC extensions seem to disable abs () overloads in <cmath>.
#ifdef _MSC_EXTENSIONS
#define BOOST_UBLAS_CMATH_BAD_STD
#endif

// We must disable element proxies as they require template partial specialisation
#define BOOST_UBLAS_NO_ELEMENT_PROXIES

// This seems to be a problem in boost.config, but won't be fixed.
#ifdef __SGI_STL_PORT
#define BOOST_MSVC_STD_ITERATOR
#endif

#endif

// Version 7.1
#if BOOST_MSVC == 1310

// Seeming still some problems with function templates
#define BOOST_UBLAS_MSVC71_FUNCTION_TEMPLATE_ORDERING
// One of these workarounds is needed for MSVC 7.1 AFAIK
// (thanks to John Maddock and Martin Lauer).
// The second workaround looks like BOOST_UBLAS_QUALIFIED_TYPENAME.
#if !(defined(BOOST_UBLAS_NO_NESTED_CLASS_RELATION) || defined(BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION))
#define BOOST_UBLAS_NO_NESTED_CLASS_RELATION
#endif

#endif

#endif


// GNU Compiler Collection
#if defined (__GNUC__) && ! defined (BOOST_STRICT_CONFIG)

#if __GNUC__ >= 3 && __GNUC_MINOR__ >= 4
// By ABI definition see GCC bug id 9982
#define BOOST_UBLAS_USEFUL_ARRAY_PLACEMENT_NEW
#endif

#if __GNUC__ <= 2 && __GNUC_MINOR__ <= 95
#define BOOST_UBLAS_NO_PROXY_SHORTCUTS
#define BOOST_UBLAS_NO_MEMBER_FRIENDS
// Cannot overload basic_stream
#define BOOST_UBLAS_USE_STREAM
// Cannot access private members from member class
#define BOOST_UBLAS_NESTED_CLASS_DR45
#endif

#endif


// Intel Compiler
#if defined (BOOST_INTEL) && ! defined (BOOST_STRICT_CONFIG)

#if (BOOST_INTEL_LINUX >= 700)
// By inspection of compiler results
#define BOOST_UBLAS_USEFUL_ARRAY_PLACEMENT_NEW
#endif

#if (BOOST_INTEL < 800)
// Base traits templates syntax untested
#define BOOST_UBLAS_NO_ITERATOR_BASE_TRAITS
#endif

// Define swap for index_pair and triple.
#if (BOOST_INTEL <= 800)
namespace boost { namespace numeric { namespace ublas {
    template<class C, class IC>
    class indexed_iterator;

    template<class V>
    class index_pair;
    template<class M>
    class index_triple;
}}}

namespace std {
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
    // iter_swap also needed for ICC on Itanium?
    template<class C, class IC>
    inline
    void iter_swap (boost::numeric::ublas::indexed_iterator<C, IC> it1,
                    boost::numeric::ublas::indexed_iterator<C, IC> it2) {
        swap (*it1, *it2);
    }
}
#endif

#endif


// Comeau compiler - thanks to Kresimir Fresl
#if defined (__COMO__) && ! defined (BOOST_STRICT_CONFIG)

// Missing std::abs overloads for float types in <cmath> are in <cstdlib>
#if defined(__LIBCOMO__) && (__LIBCOMO_VERSION__ <= 31)
#include <cstdlib>
#endif

#endif


//  HP aCC C++ compiler
#if defined (__HP_aCC) && ! defined (BOOST_STRICT_CONFIG)
// No version specific configuration
// #if (__HP_aCC <= 53800)
#endif


//  SGI MIPSpro C++ compiler
#if defined (__sgi) && ! defined (BOOST_STRICT_CONFIG)

// Missing std::abs overloads for float types in <cmath> are in <cstdlib>
// This should should be library version specific.
#include <cstdlib>

#if __COMPILER_VERSION >=650
// By inspection of compiler results - thanks to Peter Schmitteckert
#define BOOST_UBLAS_USEFUL_ARRAY_PLACEMENT_NEW
#endif

#endif


// Metrowerks Codewarrior
#if defined (__MWERKS__) && ! defined (BOOST_STRICT_CONFIG)

// 8.x
#if __MWERKS__ <= 0x3003
#define BOOST_UBLAS_NO_MEMBER_FRIENDS
// Problems with basic_stream imbue
#define BOOST_UBLAS_USE_STREAM
// Base traits templates syntax untested
#define BOOST_UBLAS_NO_ITERATOR_BASE_TRAITS
// Problems type check debug functor syntax
#define BOOST_UBLAS_TYPE_CHECK 0
#endif

#endif



// Default configuration without compiler problems
#ifndef BOOST_UBLAS_TYPENAME
#define BOOST_UBLAS_TYPENAME typename
#endif
#ifndef BOOST_UBLAS_USING
#define BOOST_UBLAS_USING using
#endif


// Enable performance options in RELEASE mode
#ifdef NDEBUG

#ifndef BOOST_UBLAS_INLINE
#define BOOST_UBLAS_INLINE inline
#endif

// Do not check sizes!
#define BOOST_UBLAS_USE_FAST_SAME

// NO runtime error checks with BOOST_UBLAS_CHECK macro
#ifndef BOOST_UBLAS_CHECK_ENABLE
#define BOOST_UBLAS_CHECK_ENABLE 0
#endif

// NO type compatibility numeric checks
#ifndef BOOST_UBLAS_TYPE_CHECK
#define BOOST_UBLAS_TYPE_CHECK 0
#endif


// Disable performance options in DEBUG mode
#else

#ifndef BOOST_UBLAS_INLINE
#define BOOST_UBLAS_INLINE
#endif

// Enable runtime error checks with BOOST_UBLAS_CHECK macro. Check bounds etc
#ifndef BOOST_UBLAS_CHECK_ENABLE
#define BOOST_UBLAS_CHECK_ENABLE 1
#endif

// Type compatibiltity numeric checks
#ifndef BOOST_UBLAS_TYPE_CHECK
#define BOOST_UBLAS_TYPE_CHECK 1
#endif

#endif


/*
 * Type compatibility checks
 *  Control type compatibility numeric runtime checks for non dense matrices.
 *  Require additional storage and complexity
 */
#if BOOST_UBLAS_TYPE_CHECK
template <class Dummy>
struct disable_type_check
{
    static bool value;
};
template <class Dummy>
bool disable_type_check<Dummy>::value = false;
#endif
#ifndef BOOST_UBLAS_TYPE_CHECK_EPSILON
#define BOOST_UBLAS_TYPE_CHECK_EPSILON (type_traits<real_type>::sqrt (std::numeric_limits<real_type>::epsilon ()))
#endif
#ifndef BOOST_UBLAS_TYPE_CHECK_MIN
#define BOOST_UBLAS_TYPE_CHECK_MIN (type_traits<real_type>::sqrt ( (std::numeric_limits<real_type>::min) ()))
#endif


/*
 * General Configuration
 */

// Operator () is alreadly heavily over used. Use project member function for proxies
//#define BOOST_UBLAS_NO_PROXY_SHORTCUTS

// In order to simplify debugging is is possible to simplify expression template
// so they are restricted to a single operation
// #define BOOST_UBLAS_SIMPLE_ET_DEBUG

// Select stream types defined for IO
#if !defined(BOOST_UBLAS_USE_STREAM) && !defined(BOOST_UBLAS_USE_BASIC_STREAM)
#define BOOST_UBLAS_USE_BASIC_STREAM
#endif

// Use invariant hoisting.
// #define BOOST_UBLAS_USE_INVARIANT_HOISTING

// Use Duff's device in element access loops
// #define BOOST_UBLAS_USE_DUFF_DEVICE

// Choose evaluation method for dense vectors and matrices
#if !(defined(BOOST_UBLAS_USE_INDEXING) || defined(BOOST_UBLAS_USE_ITERATING))
#define BOOST_UBLAS_USE_INDEXING
#endif
// #define BOOST_UBLAS_USE_ITERATING
// #define BOOST_UBLAS_ITERATOR_THRESHOLD 0
// #define BOOST_UBLAS_ITERATOR_THRESHOLD ((std::numeric_limits<std::ptrdiff_t>::max) ())

// Use indexed iterators.
// #define BOOST_UBLAS_USE_INDEXED_ITERATOR

// Alignment of bounded arrays. align(16) possibly useful for ICC
#ifndef BOOST_UBLAS_BOUNDED_ARRAY_ALIGN
#define BOOST_UBLAS_BOUNDED_ARRAY_ALIGN
#endif

// Enable assignment of non conformant proxies
#define BOOST_UBLAS_NON_CONFORMANT_PROXIES

// Enable different sparse element proxies
#ifndef BOOST_UBLAS_NO_ELEMENT_PROXIES
// Sparse proxies prevent reference invalidation problems in expressions such as:
// a [1] = a [0] = 1        Thanks to Marc Duflot for spotting this.
// #define BOOST_UBLAS_STRICT_MAP_ARRAY
#define BOOST_UBLAS_STRICT_VECTOR_SPARSE
#define BOOST_UBLAS_STRICT_MATRIX_SPARSE
// Hermitian matrices use element proxies to allow assignment to conjugate triangle
#define BOOST_UBLAS_STRICT_HERMITIAN
#endif



// Define to enable compile time const propagation for reference, proxy and closure types
#define BOOST_UBLAS_CT_REFERENCE_BASE_TYPEDEFS
#define BOOST_UBLAS_CT_PROXY_BASE_TYPEDEFS
#define BOOST_UBLAS_CT_PROXY_CLOSURE_TYPEDEFS

// Define to configure special settings for reference returning members
// #define BOOST_UBLAS_REFERENCE_CONST_MEMBER
// #define BOOST_UBLAS_PROXY_CONST_MEMBER


// Include declerations and functions
#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/ublas/definitions.hpp>


#endif
