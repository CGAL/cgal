//  (C) Copyright Gennadiy Rozental 2001-2003.
//  Use, modification, and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : support for automated test cases registration mechanism
//                for simple function based test cases
// ***************************************************************************

#ifndef BOOST_AUTO_UNIT_TEST_HPP
#define BOOST_AUTO_UNIT_TEST_HPP

// Boost.Test
#include <boost/test/unit_test.hpp>

// ************************************************************************** //
// **************           auto_unit_test_registrar           ************** //
// ************************************************************************** //

namespace boost {
namespace unit_test_framework {
namespace detail {

inline boost::unit_test_framework::test_suite*
auto_unit_test_suite()
{
    static boost::unit_test_framework::test_suite* inst = BOOST_TEST_SUITE( "Auto Unit Test" );

    return inst;
}

struct auto_unit_test_registrar
{
    // Constructor
    explicit auto_unit_test_registrar( test_case* tc ) { auto_unit_test_suite()->add( tc ); }
};

} // detail

} // unit_test_framework

} // namespace boost

// ************************************************************************** //
// **************             BOOST_AUTO_UNIT_TEST             ************** //
// ************************************************************************** //

#define BOOST_AUTO_UNIT_TEST( func_name )                           \
static void func_name();                                            \
static boost::unit_test_framework::detail::auto_unit_test_registrar \
    BOOST_JOIN( test_registrar, __LINE__)                           \
        ( BOOST_TEST_CASE( func_name ) );                           \
static void func_name()                                             \
/**/

// ************************************************************************** //
// **************             BOOST_AUTO_UNIT_TEST             ************** //
// ************************************************************************** //

#ifdef BOOST_AUTO_TEST_MAIN
boost::unit_test_framework::test_suite*
init_unit_test_suite( int /* argc */, char* /* argv */ [] ) {
    return boost::unit_test_framework::detail::auto_unit_test_suite();
}
#endif

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/05/23 10:51:33  spion
//  Initial revision
//
//  Revision 1.7  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//

// ***************************************************************************

#endif // BOOST_AUTO_UNIT_TEST_HPP
