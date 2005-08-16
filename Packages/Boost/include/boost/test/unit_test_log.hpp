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
//  Description : defines singleton class unit_test_log and all manipulators.
//  unit_test_log has output stream like interface. It's implementation is
//  completely hidden with pimple idiom
// ***************************************************************************

#ifndef BOOST_TEST_UNIT_TEST_LOG_HPP_071894GER
#define BOOST_TEST_UNIT_TEST_LOG_HPP_071894GER

// Boost.Test
#include <boost/test/test_observer.hpp>

#include <boost/test/detail/global_typedef.hpp>
#include <boost/test/detail/log_level.hpp>
#include <boost/test/detail/fwd_decl.hpp>

#include <boost/test/utils/trivial_singleton.hpp>

// Boost
#include <boost/utility.hpp>

// STL
#include <iosfwd>   // for std::ostream&

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

namespace boost {

namespace unit_test {

// ************************************************************************** //
// **************                log manipulators              ************** //
// ************************************************************************** //

namespace log {

struct begin {};

struct end {};

struct line {
    explicit    line( std::size_t ln ) : m_line_num( ln ) {}

    std::size_t m_line_num;
};

struct file {
    explicit    file( const_string fn ) : m_file_name( fn ) {}

    const_string m_file_name;
};

struct checkpoint {
    explicit    checkpoint( const_string message ) : m_message( message ) {}

    const_string m_message;
};

} // namespace log

// ************************************************************************** //
// **************             entry_value_collector            ************** //
// ************************************************************************** //

namespace ut_detail {

class entry_value_collector {
public:
    // Constructors
    entry_value_collector() : m_last( true ) {}
    entry_value_collector( entry_value_collector& rhs ) : m_last( true ) { rhs.m_last = false; }
    ~entry_value_collector();

    // collection interface
    entry_value_collector operator<<( const_string );
    entry_value_collector operator<<( log::checkpoint const& );

private:
    // Data members
    bool    m_last;
};

} // namespace ut_detail

// ************************************************************************** //
// **************                 unit_test_log                ************** //
// ************************************************************************** //

class unit_test_log_t : public test_observer, public singleton<unit_test_log_t> {
public:
    // test_observer interface implementation
    void                test_start( counter_t test_cases_amount );
    void                test_finish();
    void                test_aborted();

    void                test_unit_start( test_unit const& );
    void                test_unit_finish( test_unit const&, unsigned long elapsed );
    void                test_unit_skipped( test_unit const& );
    void                test_unit_aborted( test_unit const& );

    void                assertion_result( bool passed );
    void                exception_caught( execution_exception const& );

    // log configuration methods
    void                set_stream( std::ostream& );
    void                set_threshold_level( log_level );
    void                set_format( output_format );
    void                set_formatter( unit_test_log_formatter* );

    // entry logging
    unit_test_log_t&    operator<<( log::begin const& );        // begin entry 
    unit_test_log_t&    operator<<( log::end const& );          // end entry
    unit_test_log_t&    operator<<( log::file const& );         // set entry file name
    unit_test_log_t&    operator<<( log::line const& );         // set entry line number
    unit_test_log_t&    operator<<( log::checkpoint const& );   // set checkpoint
    unit_test_log_t&    operator<<( log_level );                // set entry level
    unit_test_log_t&    operator<<( const_string );             // log entry value

    ut_detail::entry_value_collector operator()( log_level );   // initiate entry collection

private:
    BOOST_TEST_SINGLETON_CONS( unit_test_log_t );
}; // unit_test_log_t

BOOST_TEST_SINGLETON_INST( unit_test_log )

// helper macros
#define BOOST_UT_LOG_ENTRY                                             \
    (boost::unit_test::unit_test_log << boost::unit_test::log::begin() \
        << boost::unit_test::log::file( BOOST_TEST_L( __FILE__ ) )     \
        << boost::unit_test::log::line( __LINE__ ))                    \
/**/

} // namespace unit_test

} // namespace boost

//____________________________________________________________________________//

#include <boost/test/detail/enable_warnings.hpp>

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.3  2005/08/16 11:24:12  spion
//  Import of Boost v. 1.33.0
//
//  Revision 1.30  2005/02/20 08:27:06  rogeeff
//  This a major update for Boost.Test framework. See release docs for complete list of fixes/updates
//
//  Revision 1.29  2005/02/02 12:08:14  rogeeff
//  namespace log added for log manipulators
//
//  Revision 1.28  2005/02/01 06:40:06  rogeeff
//  copyright update
//  old log entries removed
//  minor stilistic changes
//  depricated tools removed
//
//  Revision 1.27  2005/01/30 03:26:29  rogeeff
//  return an ability for explicit end()
//
//  Revision 1.26  2005/01/21 07:30:24  rogeeff
//  to log testing time log formatter interfaces changed
//
//  Revision 1.25  2005/01/18 08:26:12  rogeeff
//  unit_test_log rework:
//     eliminated need for ::instance()
//     eliminated need for << end and ...END macro
//     straitend interface between log and formatters
//     change compiler like formatter name
//     minimized unit_test_log interface and reworked to use explicit calls
//
// ***************************************************************************

#endif // BOOST_TEST_UNIT_TEST_LOG_HPP_071894GER

