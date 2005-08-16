//  (C) Copyright Gennadiy Rozental 2001-2005.
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

#ifndef BOOST_TEST_AUTO_UNIT_TEST_HPP_071894GER
#define BOOST_TEST_AUTO_UNIT_TEST_HPP_071894GER

// Boost.Test
#include <boost/test/unit_test.hpp>

#include <boost/test/detail/suppress_warnings.hpp>

// STL
#include <list>

//____________________________________________________________________________//

// ************************************************************************** //
// **************           auto_test_unit_registrar           ************** //
// ************************************************************************** //

namespace boost {
namespace unit_test {

struct auto_unit_test_suite_t : test_suite {
    auto_unit_test_suite_t() 
    : test_suite( "Master Test Suite" )
    , argc( 0 )
    , argv( 0 )
    {}
    
    // Data members    
    int      argc;
    char**   argv;
};

//____________________________________________________________________________//

inline auto_unit_test_suite_t*
auto_unit_test_suite()
{
    static auto_unit_test_suite_t* inst = new auto_unit_test_suite_t;

    return inst;
}

//____________________________________________________________________________//

namespace ut_detail {

struct auto_test_unit_registrar
{
    // Constructor
    explicit    auto_test_unit_registrar( test_case* tc, counter_t exp_fail )
    {
        curr_ts_store().back()->add( tc, exp_fail );
    }
    explicit    auto_test_unit_registrar( test_suite* ts )
    {
        curr_ts_store().back()->add( ts );

        curr_ts_store().push_back( ts );
    }
    explicit    auto_test_unit_registrar( test_unit_generator const& tc_gen )
    {
        curr_ts_store().back()->add( tc_gen );
    }
    explicit    auto_test_unit_registrar( int )
    {
        if( curr_ts_store().size() > 1 )
            curr_ts_store().pop_back();
        // else report error
    }

private:
    static std::list<test_suite*>& curr_ts_store()
    {
        static std::list<test_suite*> inst( 1, auto_unit_test_suite() );
        return inst;
    }
};

//____________________________________________________________________________//

template<typename T>
struct auto_tc_exp_fail {
    enum { value = 0 };
};

} // namespace ut_detail

} // namespace unit_test
} // namespace boost

#define BOOST_AUTO_TC_REGISTRAR( test_name )    \
    static boost::unit_test::ut_detail::auto_test_unit_registrar BOOST_JOIN( test_name, _registrar )
#define BOOST_AUTO_TC_INVOKER( test_name )      BOOST_JOIN( test_name, _invoker )
#define BOOST_AUTO_TC_UNIQUE_ID( test_name )    BOOST_JOIN( test_name, _id )

// ************************************************************************** //
// **************             BOOST_AUTO_TEST_SUITE            ************** //
// ************************************************************************** //

#define BOOST_AUTO_TEST_SUITE( suite_name )                             \
BOOST_AUTO_TC_REGISTRAR( suite_name )( BOOST_TEST_SUITE(                \
    BOOST_STRINGIZE( suite_name ) ) )                                   \
/**/

// ************************************************************************** //
// **************           BOOST_AUTO_TEST_SUITE_END          ************** //
// ************************************************************************** //

#define BOOST_AUTO_TEST_SUITE_END()                                     \
BOOST_AUTO_TC_REGISTRAR( BOOST_JOIN( end_suite, __LINE__ ) )( 1 )       \
/**/

// ************************************************************************** //
// **************             BOOST_AUTO_TEST_CASE             ************** //
// ************************************************************************** //

#define BOOST_AUTO_TEST_CASE( test_name )                               \
struct BOOST_AUTO_TC_UNIQUE_ID( test_name ) {};                         \
                                                                        \
static void test_name();                                                \
BOOST_AUTO_TC_REGISTRAR( test_name )( BOOST_TEST_CASE( test_name ),     \
    boost::unit_test::ut_detail::auto_tc_exp_fail<                      \
        BOOST_AUTO_TC_UNIQUE_ID( test_name )>::value );                 \
static void test_name()                                                 \
/**/

// ************************************************************************** //
// **************    BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES    ************** //
// ************************************************************************** //

#define BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES( test_name, n )          \
struct BOOST_AUTO_TC_UNIQUE_ID( test_name );                            \
namespace boost { namespace unit_test { namespace ut_detail {           \
                                                                        \
template<>                                                              \
struct auto_tc_exp_fail<BOOST_AUTO_TC_UNIQUE_ID( test_name ) > {        \
    enum { value = n };                                                 \
};                                                                      \
                                                                        \
}}}                                                                     \
/**/

// ************************************************************************** //
// **************            BOOST_FIXTURE_TEST_CASE           ************** //
// ************************************************************************** //

#define BOOST_FIXTURE_TEST_CASE( test_name, F )                         \
struct test_name : public F { void test_method(); };                    \
                                                                        \
void BOOST_AUTO_TC_INVOKER( test_name )()                               \
{                                                                       \
    test_name t;                                                        \
    t.test_method();                                                    \
}                                                                       \
                                                                        \
struct BOOST_AUTO_TC_UNIQUE_ID( test_name ) {};                         \
                                                                        \
BOOST_AUTO_TC_REGISTRAR( test_name )(                                   \
    boost::unit_test::make_test_case(                                   \
        &BOOST_AUTO_TC_INVOKER( test_name ), #test_name ),              \
    boost::unit_test::ut_detail::auto_tc_exp_fail<                      \
        BOOST_AUTO_TC_UNIQUE_ID( test_name )>::value  );                \
                                                                        \
void test_name::test_method()                                           \
/**/

// ************************************************************************** //
// **************        BOOST_AUTO_TEST_CASE_TEMPLATE         ************** //
// ************************************************************************** //

#define BOOST_AUTO_TEST_CASE_TEMPLATE( test_name, type_name, TL )       \
template<typename type_name>                                            \
void test_name( boost::type<type_name>* );                              \
                                                                        \
struct BOOST_AUTO_TC_INVOKER( test_name ) {                             \
    template<typename TestType>                                         \
    static void run( boost::type<TestType>* frwrd = 0 )                 \
    {                                                                   \
       test_name( frwrd );                                              \
    }                                                                   \
};                                                                      \
                                                                        \
BOOST_AUTO_TC_REGISTRAR( test_nase )(                                   \
    boost::unit_test::ut_detail::template_test_case_gen<                \
        BOOST_AUTO_TC_INVOKER( test_name ),TL >(                        \
          BOOST_STRINGIZE( test_name ) ) );                             \
                                                                        \
template<typename type_name>                                            \
void test_name( boost::type<type_name>* )                               \
/**/

// ************************************************************************** //
// **************             BOOST_AUTO_TEST_MAIN             ************** //
// ************************************************************************** //

#ifdef BOOST_AUTO_TEST_MAIN
boost::unit_test::test_suite*
init_unit_test_suite( int argc, char* argv[] ) {
    boost::unit_test::auto_unit_test_suite_t* master_test_suite = boost::unit_test::auto_unit_test_suite();

    boost::unit_test::const_string new_name = boost::unit_test::const_string( BOOST_AUTO_TEST_MAIN );

    if( !new_name.is_empty() )
        boost::unit_test::assign_op( master_test_suite->p_name.value, new_name, 0 );

    master_test_suite->argc = argc;
    master_test_suite->argv = argv;

    return master_test_suite;
}
#endif

//____________________________________________________________________________//

// deprecated
#define BOOST_AUTO_UNIT_TEST( f ) BOOST_AUTO_TEST_CASE( f )

//____________________________________________________________________________//

#include <boost/test/detail/enable_warnings.hpp>

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.3  2005/08/16 11:24:10  spion
//  Import of Boost v. 1.33.0
//
//  Revision 1.17  2005/05/08 08:55:00  rogeeff
//  typos and missing descriptions fixed
//
//  Revision 1.16  2005/05/03 03:38:35  rogeeff
//  bug in fixture test cases fixed
//
//  Revision 1.15  2005/04/18 04:54:36  rogeeff
//  Major rework in auto unit test facilities\n1. auto test suite ability introduced\n2.fixtures abilities introduced\n3. Expected failures support\n4. Master test suite renaming support
//
//  Revision 1.14  2005/03/22 06:56:13  rogeeff
//  provided access to argc/argv in auto facilities
//
//  Revision 1.13  2005/02/20 08:27:05  rogeeff
//  This a major update for Boost.Test framework. See release docs for complete list of fixes/updates
//
//  Revision 1.12  2005/02/01 06:40:06  rogeeff
//  copyright update
//  old log entries removed
//  minor stilistic changes
//  depricated tools removed
//
// ***************************************************************************

#endif // BOOST_TEST_AUTO_UNIT_TEST_HPP_071894GER
