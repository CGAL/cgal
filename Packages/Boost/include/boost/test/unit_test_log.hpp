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
//  Description : defines singleton class unit_test_log and all manipulators.
//  unit_test_log has output stream like interface. It's implementation is
//  completely hidden with pimple idiom
// ***************************************************************************

#ifndef BOOST_UNIT_TEST_LOG_HPP
#define BOOST_UNIT_TEST_LOG_HPP

// Boost.Test
#include <boost/test/detail/unit_test_config.hpp>

// BOOST
#include <boost/utility.hpp>

// STL
#include <iosfwd>   // for std::ostream&
#include <string>   // for std::string&; in fact need only forward declaration

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4512) // assignment operator could not be generated
#endif

namespace boost {

namespace unit_test_framework {

//  each log level includes all subsequent higher loging levels
enum            log_level {
    log_successful_tests     = 0,
    log_test_suites          = 1,
    log_messages             = 2,
    log_warnings             = 3,
    log_all_errors           = 4, // reported by unit test macros
    log_cpp_exception_errors = 5, // uncaught C++ exceptions
    log_system_errors        = 6, // including timeouts, signals, traps
    log_fatal_errors         = 7, // including unit test macros or
                                     // fatal system errors
    log_progress_only        = 8, // only unit test progress to be reported
    log_nothing              = 9
};

// ************************************************************************** //
// **************                log_entry_data                ************** //
// ************************************************************************** //

struct log_entry_data
{
    std::string     m_file;
    std::size_t     m_line;
    log_level       m_level;

    void clear()
    {
        m_file    = std::string();
        m_line    = 0;
        m_level   = log_nothing;
    }
};

// ************************************************************************** //
// **************                checkpoint_data               ************** //
// ************************************************************************** //

struct log_checkpoint_data
{
    std::string     m_file;
    std::size_t     m_line;
    std::string     m_message;

    void clear()
    {
        m_file    = std::string();
        m_line    = 0;
        m_message = std::string();
    }
};

// ************************************************************************** //
// **************                log manipulators              ************** //
// ************************************************************************** //

struct begin {
};

struct end {
};

struct level {
    explicit    level( log_level l_ ) : m_level( l_ ) {}

    log_level m_level;
};

struct line {
    explicit    line( std::size_t ln_ ) : m_line_num( ln_ ) {}

    std::size_t m_line_num;
};

struct file {
    explicit    file( c_string_literal fn_ ) : m_file_name( fn_ ) {}

    c_string_literal m_file_name;
};

struct checkpoint {
    explicit    checkpoint( std::string const& message_ ) : m_message( message_ ) {}

    std::string const& m_message;
};

struct log_exception {
    explicit    log_exception( c_string_literal what_ ) : m_what( what_ ) {}

    c_string_literal     m_what;
};

struct log_progress {
};

// ************************************************************************** //
// **************                 unit_test_log                ************** //
// ************************************************************************** //

class test_case;
class unit_test_log_formatter;

class unit_test_log : private boost::noncopyable { //!! Singleton
public:
    // Destructor
    ~unit_test_log();

    // instance access method;
    static unit_test_log& instance();


    void            start( bool print_build_info_ = false );
    void            header( unit_test_counter test_cases_amount_ );
    void            finish( unit_test_counter test_cases_amount_ );

    // log configuration methods
    void            set_log_stream( std::ostream& str_ );
    void            set_log_threshold_level( log_level lev_ );
    void            set_log_threshold_level_by_name( std::string const& lev_ );
    void            set_log_format( std::string const& of );
    void            set_log_formatter( unit_test_log_formatter* the_formatter );
    void            clear_checkpoint();

    // test case scope tracking
    void            track_test_case_scope( test_case const& tc, bool in_out );

    // entry configuration methods
    unit_test_log&  operator<<( begin const& );         // begin entry 
    unit_test_log&  operator<<( end const& );           // end entry
    unit_test_log&  operator<<( file const& );          // set file name
    unit_test_log&  operator<<( line const& );          // set line number
    unit_test_log&  operator<<( level const& );         // set entry level
    unit_test_log&  operator<<( checkpoint const& );    // set checkpoint

    // print value_ methods
    unit_test_log&  operator<<( log_progress const& );
    unit_test_log&  operator<<( log_exception const& );
    unit_test_log&  operator<<( c_string_literal value_ );
    unit_test_log&  operator<<( std::string const& value_ );

private:
    // formatters interface
    friend class unit_test_log_formatter;
    log_entry_data      const& entry_data() const;
    log_checkpoint_data const& checkpoint_data() const;

private:
    // Constructor
    unit_test_log();

    struct          Impl;
    Impl*           m_pimpl;
}; // unit_test_log

// helper macros
#define BOOST_UT_LOG_BEGIN( file_name, line_num, loglevel )                             \
    boost::unit_test_framework::unit_test_log::instance()                               \
                                     << boost::unit_test_framework::begin()             \
                                     << boost::unit_test_framework::level( loglevel )   \
                                     << boost::unit_test_framework::file( file_name )   \
                                     << boost::unit_test_framework::line( line_num ) << \
/**/
#define BOOST_UT_LOG_END             << boost::unit_test_framework::end();

// ************************************************************************** //
// **************            test_case_scope_tracker           ************** //
// ************************************************************************** //

struct test_case_scope_tracker {
    explicit            test_case_scope_tracker( test_case const& tc ) 
    : m_tc( tc )                                    { unit_test_log::instance().track_test_case_scope( m_tc, true ); }
                        ~test_case_scope_tracker()  { unit_test_log::instance().track_test_case_scope( m_tc, false ); }

private:
    test_case const&    m_tc;
};

} // namespace unit_test_framework

} // namespace boost

#ifdef BOOST_MSVC
# pragma warning(default: 4512) // assignment operator could not be generated
# pragma warning(pop)
#endif

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/05/23 10:51:36  spion
//  Initial revision
//
//  Revision 1.20  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//

// ***************************************************************************

#endif // BOOST_UNIT_TEST_LOG_HPP

