//  (C) Copyright Gennadiy Rozental 2003.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : 
// ***************************************************************************

#ifndef BOOST_UNIT_TEST_LOG_FORMATTER_HPP_071894GER
#define BOOST_UNIT_TEST_LOG_FORMATTER_HPP_071894GER

// Boost.Test
#include <boost/test/detail/unit_test_config.hpp>
#include <boost/test/unit_test_log.hpp>

// BOOST

// STL
#include <iosfwd>
#include <string> // need only forward decl

#include <boost/test/detail/suppress_warnings.hpp>

namespace boost {

namespace unit_test {

// ************************************************************************** //
// **************            unit_test_log_formatter           ************** //
// ************************************************************************** //

class unit_test_log_formatter {
public:
    enum log_entry_types { BOOST_UTL_ET_INFO, 
                           BOOST_UTL_ET_MESSAGE,
                           BOOST_UTL_ET_WARNING,
                           BOOST_UTL_ET_ERROR,
                           BOOST_UTL_ET_FATAL_ERROR };

    // Constructor
    explicit unit_test_log_formatter( unit_test_log const& log )
    : m_log( log ) {}

    // Destructor
    virtual             ~unit_test_log_formatter() {}

    // Formatter interface
    virtual void        start_log( std::ostream& output, bool log_build_info ) = 0;
    virtual void        log_header( std::ostream& output, unit_test_counter test_cases_amount ) = 0;
    virtual void        finish_log( std::ostream& output ) = 0;

    virtual void        track_test_case_scope( std::ostream& output, test_case const& tc, bool in_out ) = 0;
    virtual void        log_exception( std::ostream& output, const_string test_case_name, const_string explanation ) = 0;

    virtual void        begin_log_entry( std::ostream& output, log_entry_types let ) = 0;
    virtual void        log_entry_value( std::ostream& output, const_string value ) = 0;
    virtual void        end_log_entry( std::ostream& output ) = 0;

protected:
    // Implementation interface
    log_entry_data      const& entry_data() const       { return m_log.entry_data(); }
    log_checkpoint_data const& checkpoint_data() const  { return m_log.checkpoint_data(); }

private:
    // Data members
    unit_test_log const& m_log;
};

} // namespace unit_test

} // namespace boost

#include <boost/test/detail/enable_warnings.hpp>

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.2  2004/11/20 10:52:16  spion
//  Import of Boost v. 1.32.0
//
//  Revision 1.7  2004/07/19 12:16:23  rogeeff
//  guard rename
//  warning suppressed
//
//  Revision 1.6  2004/05/21 06:19:35  rogeeff
//  licence update
//
//  Revision 1.5  2004/05/11 11:00:51  rogeeff
//  basic_cstring introduced and used everywhere
//  class properties reworked
//
//  Revision 1.4  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//
// ***************************************************************************

#endif // BOOST_UNIT_TEST_LOG_FORMATTER_HPP_071894GER

