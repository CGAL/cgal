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
//  Description : defines specific version of execution monitor used to run unit 
//  test cases. Translates executioin exception into error level
// ***************************************************************************

#ifndef BOOST_UNIT_TEST_MONITOR_HPP_071894GER
#define BOOST_UNIT_TEST_MONITOR_HPP_071894GER

#include <boost/test/execution_monitor.hpp>

#include <boost/test/detail/suppress_warnings.hpp>

namespace boost {

namespace unit_test {

class test_case;

namespace ut_detail {

// ************************************************************************** //
// **************               unit_test_monitor              ************** //
// ************************************************************************** //

class unit_test_monitor : public execution_monitor {
    typedef void (test_case::*function_to_monitor)();
public:
    enum error_level { 
        test_fail               =  1,
        test_ok                 =  0,
        constructor_error       = -1, 
        unexpected_exception    = -2, 
        os_exception            = -3, 
        os_timeout              = -4, 
        fatal_error             = -5,  // includes both system and user
        destructor_error        = -6
    };

    static bool         is_critical_error( error_level e_ ) { return e_ <= fatal_error; }

    // management method; same for all monitors
    static void         catch_system_errors( bool yes_no = true ) { s_catch_system_errors = yes_no; }

    // monitor method
    error_level         execute_and_translate( test_case* target_test_case_, function_to_monitor f_, int timeout_ );

    // execution monitor hook implementation
    virtual int         function();

private:
    // Data members
    function_to_monitor m_test_case_method;
    test_case*          m_test_case;
    static bool         s_catch_system_errors;
}; // unit_test_monitor

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
//  Revision 1.17  2004/07/19 12:24:01  rogeeff
//  guard rename
//  suppress warnings reworked
//
//  Revision 1.16  2004/06/07 07:33:49  rogeeff
//  detail namespace renamed
//
//  Revision 1.15  2004/05/21 06:19:35  rogeeff
//  licence update
//
//  Revision 1.14  2004/05/11 11:00:53  rogeeff
//  basic_cstring introduced and used everywhere
//  class properties reworked
//
//  Revision 1.13  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//
// ***************************************************************************

#endif // BOOST_UNIT_TEST_MONITOR_HPP_071894GER
