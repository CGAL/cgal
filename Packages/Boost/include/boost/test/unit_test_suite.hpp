//  (C) Copyright Gennadiy Rozental 2005.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : defines test_unit, test_case, test_case_results, test_suite and test_tree_visitor
// ***************************************************************************

#ifndef BOOST_TEST_UNIT_TEST_SUITE_HPP_071894GER
#define BOOST_TEST_UNIT_TEST_SUITE_HPP_071894GER

// Boost.Test
#include <boost/test/detail/config.hpp>
#include <boost/test/detail/global_typedef.hpp>
#include <boost/test/utils/class_properties.hpp>
#include <boost/test/utils/callback.hpp>
#include <boost/test/detail/fwd_decl.hpp>
#include <boost/test/detail/workaround.hpp>

// Boost
#include <boost/shared_ptr.hpp>

// STL
#include <string>   // for std::string
#include <list>     // for std::list
#include <vector>   // for std::list

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

#define BOOST_TEST_CASE( function ) \
boost::unit_test::make_test_case( boost::unit_test::callback0<>(function), BOOST_TEST_STRINGIZE( function ) )
#define BOOST_CLASS_TEST_CASE( function, tc_instance ) \
boost::unit_test::make_test_case((function), BOOST_TEST_STRINGIZE( function ), tc_instance )
#define BOOST_TEST_SUITE( testsuite_name ) \
( new boost::unit_test::test_suite( testsuite_name ) )

namespace boost {

namespace unit_test {

// ************************************************************************** //
// **************                   test_unit                  ************** //
// ************************************************************************** //

class test_unit {
public:
    enum { type = tut_any };

    // Constructor
    test_unit( const_string tu_name, test_unit_type t );

    // dependencies management
    void    depends_on( test_unit* tu );
    bool    check_dependencies() const;

    // Public r/o properties
    typedef BOOST_READONLY_PROPERTY(test_unit_id,(framework_impl)) id_t;
    readonly_property<test_unit_type>   p_type;                 // type for this test unit
    readonly_property<const_string>     p_type_name;            // "case"/"suite"
    id_t                                p_id;                   // unique id for this test unit

    // Public r/w properties
    readwrite_property<std::string>     p_name;                 // name for this test unit
    readwrite_property<unsigned>        p_timeout;              // timeout for the test unit execution 
    readwrite_property<counter_t>       p_expected_failures;    // number of expected failured in this test unit

private:
    // Data members
    std::list<test_unit_id>             m_dependencies;
};

// ************************************************************************** //
// **************              test_case_generator             ************** //
// ************************************************************************** //

class test_unit_generator {
public:
    virtual test_unit*  next() const = 0;

protected:
    BOOST_TEST_PROTECTED_VIRTUAL ~test_unit_generator() {}
};

// ************************************************************************** //
// **************                   test_case                  ************** //
// ************************************************************************** //

class test_case : public test_unit {
public:
    enum { type = tut_case };

    // Constructor
    test_case( const_string tc_name, callback0<> const& test_func );

    // Access methods
    callback0<> const&  test_func() const { return m_test_func; }

private:
    friend class framework_impl;
    ~test_case() {}

    // BOOST_MSVC <= 1200 have problems with callback as property
    // Data members
    callback0<> m_test_func;
};

// ************************************************************************** //
// **************                  test_suite                  ************** //
// ************************************************************************** //

class test_suite : public test_unit {
public:
    enum { type = tut_suite };

    // Constructor
    explicit    test_suite( const_string ts_name = "Master" );

    // test case list management
    void        add( test_unit* tu, counter_t expected_failures = 0, unsigned timeout = 0 );
    void        add( test_unit_generator const& gen, unsigned timeout = 0 );

protected:
    friend void traverse_test_tree( test_suite const&, test_tree_visitor& );
    friend class framework_impl;
    virtual     ~test_suite() {}

private:
    // Data members
    std::vector<test_unit_id> m_members;
};

// ************************************************************************** //
// **************               test_tree_visitor              ************** //
// ************************************************************************** //

class test_tree_visitor {
public:
    // test tree visitor interface
    virtual void    visit( test_case const& )               {}
    virtual bool    test_suite_start( test_suite const& )   { return true; }
    virtual void    test_suite_finish( test_suite const& )  {}

protected:
    BOOST_TEST_PROTECTED_VIRTUAL ~test_tree_visitor() {}
};

// ************************************************************************** //
// **************               traverse_test_tree             ************** //
// ************************************************************************** //

void    traverse_test_tree( test_case const&, test_tree_visitor& );
void    traverse_test_tree( test_suite const&, test_tree_visitor& );
void    traverse_test_tree( test_unit_id id, test_tree_visitor& );

//____________________________________________________________________________//

inline void
traverse_test_tree( test_unit const& tu, test_tree_visitor& V )
{
    if( tu.p_type == tut_case )
        traverse_test_tree( static_cast<test_case const&>( tu ), V );
    else
        traverse_test_tree( static_cast<test_suite const&>( tu ), V );
}

//____________________________________________________________________________//

// ************************************************************************** //
// **************                test_case_counter             ************** //
// ************************************************************************** //

struct test_case_counter : test_tree_visitor {
    test_case_counter() : m_count( 0 ) {}

    void        visit( test_case const& ) { m_count++; }

    counter_t   m_count;
};

// ************************************************************************** //
// **************                  test_aborted                ************** //
// ************************************************************************** //

struct test_aborted {};

// ************************************************************************** //
// **************               object generators              ************** //
// ************************************************************************** //

namespace ut_detail {

std::string normalize_test_case_name( const_string tu_name );

template<typename UserTestCase>
struct user_tc_method_invoker {
    typedef void (UserTestCase::*test_method )();

    user_tc_method_invoker( shared_ptr<UserTestCase> inst, test_method tm )
    : m_inst( inst ), m_test_method( tm ) {}

    void operator()() { ((*m_inst).*m_test_method)(); }

    shared_ptr<UserTestCase> m_inst;
    test_method              m_test_method;
};

} // namespace ut_detail

//____________________________________________________________________________//

inline test_case*
make_test_case( callback0<> const& test_func, const_string tc_name )
{
    return new test_case( ut_detail::normalize_test_case_name( tc_name ), test_func );
}

//____________________________________________________________________________//

template<typename UserTestCase>
inline test_case*
make_test_case( void (UserTestCase::*test_method )(),
                  const_string tc_name,
                  boost::shared_ptr<UserTestCase> const& user_test_case )
{
    return new test_case( ut_detail::normalize_test_case_name( tc_name ), 
                          ut_detail::user_tc_method_invoker<UserTestCase>( user_test_case, test_method ) );
}

//____________________________________________________________________________//

} // unit_test

} // namespace boost

//____________________________________________________________________________//

#include <boost/test/detail/enable_warnings.hpp>

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.3  2005/08/16 11:24:10  spion
//  Import of Boost v. 1.33.0
//
//  Revision 1.32  2005/05/02 06:00:10  rogeeff
//  restore a parameterized user case method based testing
//
//  Revision 1.31  2005/04/18 04:55:30  rogeeff
//  test unit name made read/write
//
//  Revision 1.30  2005/03/22 06:57:29  rogeeff
//  allow to inherit test_suite
//
//  Revision 1.29  2005/02/21 10:25:54  rogeeff
//  use std::vector so we could employ random_shuffle
//
//  Revision 1.28  2005/02/20 08:27:06  rogeeff
//  This a major update for Boost.Test framework. See release docs for complete list of fixes/updates
//
//  Revision 1.27  2005/02/01 06:40:06  rogeeff
//  copyright update
//  old log entries removed
//  minor stylistic changes
//  deprecated tools removed
//
//  Revision 1.26  2005/01/30 03:22:07  rogeeff
//  interface changed to use const_string
//  use BOOST_TEST_STRINGIZE
//
//  Revision 1.25  2005/01/22 19:22:12  rogeeff
//  implementation moved into headers section to eliminate dependency of included/minimal component on src directory
//
// ***************************************************************************

#endif // BOOST_TEST_UNIT_TEST_SUITE_HPP_071894GER

