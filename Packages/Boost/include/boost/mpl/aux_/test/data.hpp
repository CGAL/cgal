
#ifndef BOOST_MPL_AUX_TEST_DATA_HPP_INCLUDED
#define BOOST_MPL_AUX_TEST_DATA_HPP_INCLUDED

// + file: boost/mpl/aux_/test/data.hpp
// + last modified: 04/may/03

// Copyright (c) 2002-03
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.
//
// See http://www.boost.org/libs/mpl for documentation.

#include "boost/noncopyable.hpp"

struct UDT {};
struct incomplete;
class abstract { virtual ~abstract() = 0; };
using boost::noncopyable;

// to do: add function types for compilers that are able to handle them
// (T ())(T (int))(void (T))(T (T)) 
//    (T[10])(T const[10])(T volatile[10])(T const volatile[10])
//    (T (*)())(T (* const)())(T (* volatile)())(T (* const volatile)())
//    (T (*)(int))(void (*)(T))(T (*)(T))

#define CTT_basic_modifiers( T ) \
    (T const)(T volatile)(T const volatile) \
    (T&)(T const&)(T volatile&)(T const volatile&) \
    (T*)(T const*)(T volatile*)(T const volatile*) \
/**/

#define CTT_basic_types() \
    (bool)(char)(int)(UDT)(incomplete)(noncopyable)(abstract) \
/**/

#endif // BOOST_MPL_AUX_TEST_DATA_HPP_INCLUDED
