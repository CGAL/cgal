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
//  Description : contains log formatters supplied by the framework definitions 
// ***************************************************************************

#ifndef BOOST_TEST_SUPPLIED_LOG_FORMATTERS_HPP
#define BOOST_TEST_SUPPLIED_LOG_FORMATTERS_HPP

// Boost.Test
#include <boost/test/unit_test_log_formatter.hpp>

// BOOST
#include <boost/config.hpp>

// STL
#include <cstddef>

namespace boost {

namespace unit_test_framework {

namespace detail {

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
    void    log_exception( std::ostream& output, std::string const& test_case_name, c_string_literal explanation );
    void    begin_log_entry( std::ostream& output, log_entry_types let );

    void    log_entry_value( std::ostream& output, std::string const& value );
    void    end_log_entry( std::ostream& output );

protected:
    virtual void    print_prefix( std::ostream& output, std::string const& file, std::size_t line );
};

// ************************************************************************** //
// **************               xml_log_formatter              ************** //
// ************************************************************************** //

class xml_log_formatter : public unit_test_log_formatter {
public:
    explicit xml_log_formatter( unit_test_log const& log );

    void    start_log( std::ostream& output, bool log_build_info );
    void    log_header( std::ostream& output, unit_test_counter test_cases_amount );
    void    finish_log( std::ostream& output );

    void    track_test_case_scope( std::ostream& output, test_case const& tc, bool in_out );
    void    log_exception( std::ostream& output, std::string const& test_case_name, c_string_literal explanation );
    void    begin_log_entry( std::ostream& output, log_entry_types let );

    void    log_entry_value( std::ostream& output, std::string const& value );
    void    end_log_entry( std::ostream& output );

private:
    void    print_indent( std::ostream& output );

    // Data members
    std::size_t      m_indent;
    c_string_literal m_curr_tag;
};

} // namespace detail

} // namespace unit_test_framework

} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/05/23 10:51:39  spion
//  Initial revision
//
//  Revision 1.4  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//

// ***************************************************************************

#endif // BOOST_TEST_SUPPLIED_LOG_FORMATTERS_HPP
