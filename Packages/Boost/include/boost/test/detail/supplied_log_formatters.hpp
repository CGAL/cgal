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
//  Description : contains log formatters supplied by the framework definitions 
// ***************************************************************************

#ifndef BOOST_TEST_SUPPLIED_LOG_FORMATTERS_HPP_071894GER
#define BOOST_TEST_SUPPLIED_LOG_FORMATTERS_HPP_071894GER

// Boost.Test
#include <boost/test/unit_test_log_formatter.hpp>
#include <boost/test/detail/xml_printer.hpp>

// BOOST
#include <boost/config.hpp>

// STL
#include <cstddef>

#include <boost/test/detail/suppress_warnings.hpp>

namespace boost {

namespace unit_test {

namespace ut_detail {

// ************************************************************************** //
// **************           msvc65_like_log_formatter          ************** //
// ************************************************************************** //

class msvc65_like_log_formatter : public unit_test_log_formatter {
public:
    explicit msvc65_like_log_formatter( unit_test_log const& log );

    void    start_log( std::ostream& output, bool log_build_info );
    void    log_header( std::ostream& output, unit_test_counter test_cases_amount );
    void    finish_log( std::ostream& output );

    void    track_test_case_scope( std::ostream& output, test_case const& tc, bool in_out );
    void    log_exception( std::ostream& output, const_string test_case_name, const_string explanation );
    void    begin_log_entry( std::ostream& output, log_entry_types let );

    void    log_entry_value( std::ostream& output, const_string value );
    void    end_log_entry( std::ostream& output );

protected:
    virtual void    print_prefix( std::ostream& output, const_string file, std::size_t line );
};

// ************************************************************************** //
// **************               xml_log_formatter              ************** //
// ************************************************************************** //

class xml_log_formatter : public unit_test_log_formatter, private xml_printer {
public:
    explicit xml_log_formatter( unit_test_log const& log );

    void    start_log( std::ostream& output, bool log_build_info );
    void    log_header( std::ostream& output, unit_test_counter test_cases_amount );
    void    finish_log( std::ostream& output );

    void    track_test_case_scope( std::ostream& output, test_case const& tc, bool in_out );
    void    log_exception( std::ostream& output, const_string test_case_name, const_string explanation );
    void    begin_log_entry( std::ostream& output, log_entry_types let );

    void    log_entry_value( std::ostream& output, const_string value );
    void    end_log_entry( std::ostream& output );

private:
    void    print_indent( std::ostream& output );

    // Data members
    std::size_t     m_indent;
    const_string    m_curr_tag;
};

} // namespace ut_detail

} // namespace unit_test

} // namespace boost

#include <boost/test/detail/enable_warnings.hpp>

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.2  2004/11/20 10:52:20  spion
//  Import of Boost v. 1.32.0
//
//  Revision 1.8  2004/07/19 12:22:49  rogeeff
//  guard rename
//  suppress warnings
//
//  Revision 1.7  2004/06/07 07:33:49  rogeeff
//  detail namespace renamed
//
//  Revision 1.6  2004/05/21 06:19:35  rogeeff
//  licence update
//
//  Revision 1.5  2004/05/11 11:00:53  rogeeff
//  basic_cstring introduced and used everywhere
//  class properties reworked
//
//  Revision 1.4  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//
// ***************************************************************************

#endif // BOOST_TEST_SUPPLIED_LOG_FORMATTERS_HPP_071894GER
