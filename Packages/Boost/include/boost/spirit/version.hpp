/*=============================================================================
    Copyright (c) 2001-2003 Hartmut Kaiser
    http://spirit.sourceforge.net/

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#if !defined(SPIRIT_VERSION_HPP)
#define SPIRIT_VERSION_HPP

///////////////////////////////////////////////////////////////////////////////
//
//  This checks, whether the used Boost library is at least V1.32.0
//
///////////////////////////////////////////////////////////////////////////////
#include <boost/version.hpp>

#if BOOST_VERSION < 103200
#error "Spirit V1.8.1 needs at least Boost V1.32.0 to compile successfully."
#endif 

///////////////////////////////////////////////////////////////////////////////
//
//  This is the version of the current Spirit distribution
//
///////////////////////////////////////////////////////////////////////////////
#define SPIRIT_VERSION 0x1801
#define SPIRIT_PIZZA_VERSION SPIRIT_FOUR_SEASONS  // :-)

#endif // defined(SPIRIT_VERSION_HPP)
