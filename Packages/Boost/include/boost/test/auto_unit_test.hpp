//  (C) Copyright Gennadiy Rozental 2001-2004.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : support for automated test cases registration mechanism
//                for simple function based test cases
// ***************************************************************************

#ifndef BOOST_AUTO_UNIT_TEST_HPP_071894GER
#define BOOST_AUTO_UNIT_TEST_HPP_071894GER

// Boost.Test
#include <boost/test/unit_test.hpp>

// ************************************************************************** //
// **************           auto_unit_test_registrar           ************** //
// ************************************************************************** //

namespace boost {
namespace unit_test {
namespace ut_detail {

inline boost::unit_test::test_suite*
auto_unit_test_suite()
{
    static boost::unit_test::test_suite* inst = BOOST_TEST_SUITE( "Auto Unit Test" );

    return inst;
}

struct auto_unit_test_registrar
{
    // Constructor
    explicit auto_unit_test_registrar( test_case* tc ) { auto_unit_test_suite()->add( tc ); }
};

} // namespace ut_detail

} // namespace unit_test

} // namespace boost

// ************************************************************************** //
// **************             BOOST_AUTO_UNIT_TEST             ************** //
// ************************************************************************** //

#define BOOST_AUTO_UNIT_TEST( func_name )                           \
static void func_name();                                            \
static boost::unit_test::ut_detail::auto_unit_test_registrar        \
    BOOST_JOIN( test_registrar, __LINE__)                           \
        ( BOOST_TEST_CASE( func_name ) );                           \
static void func_name()                                             \
/**/

// ************************************************************************** //
// **************             BOOST_AUTO_TEST_MAIN             ************** //
// ************************************************************************** //

#ifdef BOOST_AUTO_TEST_MAIN
boost::unit_test::test_suite*
init_unit_test_suite( int /* argc */, char* /* argv */ [] ) {
    return boost::unit_test::ut_detail::auto_unit_test_suite();
}
#endif

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.2  2004/11/20 10:52:12  spion
//  Import of Boost v. 1.32.0
//
//  Revision 1.11  2004/07/19 12:12:40  rogeeff
//  guard rename
//
//  Revision 1.10  2004/06/07 07:33:42  rogeeff
//  detail namespace renamed
//
//  Revision 1.9  2004/05/21 06:19:35  rogeeff
//  licence update
//
//  Revision 1.8  2004/05/11 11:00:33  rogeeff
//  basic_cstring introduced and used everywhere
//  class properties reworked
//
//  Revision 1.7  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//
// ***************************************************************************

#endif // BOOST_AUTO_UNIT_TEST_HPP_071894GER
