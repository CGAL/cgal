
//  (C) Copyright Steve Cleary, Beman Dawes, Howard Hinnant & John Maddock 2000.
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).
//
//  See http://www.boost.org/libs/type_traits for most recent version including documentation.

#ifndef BOOST_TT_CONFIG_HPP_INCLUDED
#define BOOST_TT_CONFIG_HPP_INCLUDED

#ifndef BOOST_CONFIG_HPP
#include "boost/config.hpp"
#endif

//
// whenever we have a conversion function with elipses
// it needs to be declared __cdecl to suppress compiler
// warnings from MS and Borland compilers (this *must*
// appear before we include is_same.hpp below):
#if defined(BOOST_MSVC) || (defined(__BORLANDC__) && !defined(BOOST_DISABLE_WIN32))
#   define BOOST_TT_DECL __cdecl
#else
#   define BOOST_TT_DECL /**/
#endif

# if (defined(__MWERKS__) && __MWERKS__ >= 0x3000) || (defined(BOOST_MSVC) && (BOOST_MSVC > 1301)) || defined(__EDG_VERSION__) || (defined(__GNUC__) && (__GNUC__ >= 3)) || defined(__DMC__) || defined(BOOST_NO_COMPILER_CONFIG)
#   define BOOST_TT_HAS_CONFORMING_IS_CLASS_IMPLEMENTATION
#endif

#endif // BOOST_TT_CONFIG_HPP_INCLUDED


