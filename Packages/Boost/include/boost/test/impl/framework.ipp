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
//  Description : implements framework songleton - main driver for the test
// ***************************************************************************

#ifndef BOOST_TEST_FRAMEWORK_IPP_021005GER
#define BOOST_TEST_FRAMEWORK_IPP_021005GER

// Boost.Test
#include <boost/test/framework.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/unit_test_log.hpp>
#include <boost/test/unit_test_monitor.hpp>
#include <boost/test/test_observer.hpp>
#include <boost/test/results_collector.hpp>
#include <boost/test/progress_monitor.hpp>
#include <boost/test/results_reporter.hpp>
#include <boost/test/test_tools.hpp>

#include <boost/test/detail/unit_test_parameters.hpp>
#include <boost/test/detail/global_typedef.hpp>

#include <boost/test/utils/foreach.hpp>

// Boost
#include <boost/timer.hpp>

// STL
#include <map>
#include <stdexcept>
#include <cstdlib>
#include <ctime>

#ifdef BOOST_NO_STDC_NAMESPACE
namespace std { using ::time; using ::srand; }
#endif

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

// prototype for user's test suite init function
extern boost::unit_test::test_suite* init_unit_test_suite( int argc, char* argv[] );

namespace boost {

namespace unit_test {

// ************************************************************************** //
// **************                   framework                  ************** //
// ************************************************************************** //

class framework_impl : public test_tree_visitor {
public:
    framework_impl()
    : m_master_test_suite( INV_TEST_UNIT_ID )
    , m_curr_test_case( INV_TEST_UNIT_ID )
    , m_next_test_case_id( MIN_TEST_CASE_ID )
    , m_next_test_suite_id( MIN_TEST_SUITE_ID )
    , m_test_in_progress( false )
    {}

    ~framework_impl()
    {
        BOOST_TEST_FOREACH( test_unit_store::value_type const&, tu, m_test_units ) {
            if( test_id_2_unit_type( tu.second->p_id ) == tut_suite )
                delete  (test_suite const*)tu.second;
            else
                delete  (test_case const*)tu.second;
        }
    }

    void            set_tu_id( test_unit& tu, test_unit_id id ) { tu.p_id.value = id; }

    // test_tree_visitor interface implementation
    void            visit( test_case const& tc )
    {
        if( !tc.check_dependencies() ) {
            BOOST_TEST_FOREACH( test_observer*, to, m_observers )
                to->test_unit_skipped( tc );

            return;
        }

        BOOST_TEST_FOREACH( test_observer*, to, m_observers )
            to->test_unit_start( tc );

        boost::timer tc_timer;
        test_unit_id bkup = m_curr_test_case;
        m_curr_test_case = tc.p_id;
        unit_test_monitor_t::error_level run_result = unit_test_monitor.execute_and_translate( tc );

        unsigned long elapsed = static_cast<unsigned long>( tc_timer.elapsed() * 1e6 );

        if( unit_test_monitor.is_critical_error( run_result ) ) {
            BOOST_TEST_FOREACH( test_observer*, to, m_observers )
                to->test_aborted();
        }

        BOOST_TEST_FOREACH( test_observer*, to, m_observers )
            to->test_unit_finish( tc, elapsed );

        m_curr_test_case = bkup;

        if( unit_test_monitor.is_critical_error( run_result ) )
            throw test_aborted();
    }

    bool            test_suite_start( test_suite const& ts )
    {
        if( !ts.check_dependencies() ) {
            BOOST_TEST_FOREACH( test_observer*, to, m_observers )
                to->test_unit_skipped( ts );

            return false;
        }

        BOOST_TEST_FOREACH( test_observer*, to, m_observers )
            to->test_unit_start( ts );

        return true;
    }

    void            test_suite_finish( test_suite const& ts )
    {
        BOOST_TEST_FOREACH( test_observer*, to, m_observers )
            to->test_unit_finish( ts, 0 );
    }

    //////////////////////////////////////////////////////////////////

    typedef std::map<test_unit_id,test_unit const*> test_unit_store;
    typedef std::list<test_observer*>               observer_store;

    test_unit_id    m_master_test_suite;
    test_unit_id    m_curr_test_case;
    test_unit_store m_test_units;

    test_unit_id    m_next_test_case_id;
    test_unit_id    m_next_test_suite_id;

    bool            m_test_in_progress;

    observer_store  m_observers;

};

//____________________________________________________________________________//

namespace {

framework_impl& s_frk_impl() { static framework_impl the_inst; return the_inst; }

} // local namespace

//____________________________________________________________________________//

namespace framework {

void
init( int argc, char* argv[] )
{
    runtime_config::init( &argc, argv );

    // set the log level nad format
    unit_test_log.set_threshold_level( runtime_config::log_level() );
    unit_test_log.set_format( runtime_config::log_format() );

    // set the report level nad format
    results_reporter::set_level( runtime_config::report_level() );
    results_reporter::set_format( runtime_config::report_format() );

    register_observer( results_collector );
    register_observer( unit_test_log );

    if( runtime_config::show_progress() )
        register_observer( progress_monitor );

    if( runtime_config::detect_memory_leak() > 0 )
        detect_memory_leak( runtime_config::detect_memory_leak() );

    // init master unit test suite
    test_suite const* master_suite = init_unit_test_suite( argc, argv );
    if( !master_suite )
        throw std::logic_error( "Fail to initialize test suite" );

    s_frk_impl().m_master_test_suite = master_suite->p_id;
}

//____________________________________________________________________________//

void
register_test_unit( test_case* tc )
{
    if( tc->p_id != INV_TEST_UNIT_ID )
        throw std::logic_error( "Test case already registered" );

    test_unit_id new_id = s_frk_impl().m_next_test_case_id;

    if( new_id == MAX_TEST_CASE_ID )
        throw std::logic_error( "Too many test cases" );

    typedef framework_impl::test_unit_store::value_type map_value_type;

    s_frk_impl().m_test_units.insert( map_value_type( new_id, tc ) );
    s_frk_impl().m_next_test_case_id++;

    s_frk_impl().set_tu_id( *tc, new_id );
}

//____________________________________________________________________________//

void
register_test_unit( test_suite* ts )
{
    if( ts->p_id != INV_TEST_UNIT_ID )
        throw std::logic_error( "Test suite already registered" );

    test_unit_id new_id = s_frk_impl().m_next_test_suite_id;

    if( new_id == MAX_TEST_SUITE_ID )
        throw std::logic_error( "Too many test suites" );

    typedef framework_impl::test_unit_store::value_type map_value_type;
    s_frk_impl().m_test_units.insert( map_value_type( new_id, ts ) );
    s_frk_impl().m_next_test_suite_id++;

    s_frk_impl().set_tu_id( *ts, new_id );
}

//____________________________________________________________________________//

void
register_observer( test_observer& to )
{
    s_frk_impl().m_observers.push_back( &to );
}

//____________________________________________________________________________//

void
reset_observers()
{
    s_frk_impl().m_observers.clear();
}

//____________________________________________________________________________//

test_suite const&
master_test_suite()
{
    return get<test_suite>( s_frk_impl().m_master_test_suite );
}

//____________________________________________________________________________//

test_case const&
current_test_case()
{
    return get<test_case>( s_frk_impl().m_curr_test_case );
}

//____________________________________________________________________________//

test_unit const&
get( test_unit_id id, test_unit_type t )
{
    test_unit const* res = s_frk_impl().m_test_units[id];

    if( (res->p_type & t) == 0 )
        throw std::logic_error( "Invalid test unit type" );

    return *res;
}

//____________________________________________________________________________//

void
run( test_unit_id id, bool continue_test )
{
    if( id == INV_TEST_UNIT_ID )
        id = s_frk_impl().m_master_test_suite;

    if( id == INV_TEST_UNIT_ID )
        throw std::logic_error( "Test unit is initialized" );

    test_case_counter tcc;
    traverse_test_tree( id, tcc );

    bool call_start_finish = !continue_test || !s_frk_impl().m_test_in_progress;
    bool was_in_progress = s_frk_impl().m_test_in_progress;

    s_frk_impl().m_test_in_progress = true;

    if( call_start_finish ) {
        BOOST_TEST_FOREACH( test_observer*, to, s_frk_impl().m_observers )
            to->test_start( tcc.m_count );
    }

    switch( runtime_config::random_seed() ) {
    case 0:
        break;
    case 1: {
        unsigned int seed = std::time( 0 );
        BOOST_MESSAGE( "Test cases order is shuffled using seed: " << seed );
        std::srand( seed );
        break;
    }
    default:
        BOOST_MESSAGE( "Test cases order is shuffled using seed: " << runtime_config::random_seed() );
        std::srand( runtime_config::random_seed() );
    }

    try {
        traverse_test_tree( id, s_frk_impl() );
    }
    catch( test_aborted const& ) {
        // abort already reported
    }

    if( call_start_finish ) {
        BOOST_TEST_FOREACH( test_observer*, to, s_frk_impl().m_observers )
            to->test_finish();
    }

    s_frk_impl().m_test_in_progress = was_in_progress;
}

//____________________________________________________________________________//

void
run( test_unit const* tu, bool continue_test )
{
    run( tu->p_id, continue_test );
}

//____________________________________________________________________________//

void
assertion_result( bool passed )
{
    BOOST_TEST_FOREACH( test_observer*, to, s_frk_impl().m_observers )
        to->assertion_result( passed );
}

//____________________________________________________________________________//

void
exception_caught( execution_exception const& ex )
{
    BOOST_TEST_FOREACH( test_observer*, to, s_frk_impl().m_observers )
        to->exception_caught( ex );
}

//____________________________________________________________________________//

void
test_unit_aborted()
{
    test_unit const& tu = current_test_case();

    BOOST_TEST_FOREACH( test_observer*, to, s_frk_impl().m_observers )
        to->test_unit_aborted( tu );
}

//____________________________________________________________________________//

} // namespace framework

} // namespace unit_test

} // namespace boost

//____________________________________________________________________________//

#include <boost/test/detail/enable_warnings.hpp>

// ***************************************************************************
//  Revision History :
//
//  $Log$
//  Revision 1.1  2005/08/16 11:24:13  spion
//  Initial revision
//
//  Revision 1.6  2005/05/08 08:55:09  rogeeff
//  typos and missing descriptions fixed
//
//  Revision 1.5  2005/04/05 07:23:20  rogeeff
//  restore default
//
//  Revision 1.4  2005/04/05 06:11:37  rogeeff
//  memory leak allocation point detection\nextra help with _WIN32_WINNT
//
//  Revision 1.3  2005/03/23 21:02:19  rogeeff
//  Sunpro CC 5.3 fixes
//
//  Revision 1.2  2005/02/21 10:12:18  rogeeff
//  Support for random order of test cases implemented
//
//  Revision 1.1  2005/02/20 08:27:07  rogeeff
//  This a major update for Boost.Test framework. See release docs for complete list of fixes/updates
//
// ***************************************************************************

#endif // BOOST_TEST_FRAMEWORK_IPP_021005GER
