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
//  Description : defines specific version of execution monitor used to run unit 
//  test cases. Translates executioin exception into error level
// ***************************************************************************

#ifndef BOOST_UNIT_TEST_MONITOR_HPP
#define BOOST_UNIT_TEST_MONITOR_HPP

#include <boost/test/execution_monitor.hpp>

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4511) // copy constructor could not be generated
# pragma warning(disable: 4512) // assignment operator could not be generated
#endif

namespace boost {

namespace unit_test_framework {

class test_case;

namespace detail {

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

} // namespace detail

} // namespace unit_test_framework

} // namespace boost

#ifdef BOOST_MSVC
# pragma warning(default: 4511) // copy constructor could not be generated
# pragma warning(default: 4512) // assignment operator could not be generated
# pragma warning(pop)
#endif

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/05/23 10:51:39  spion
//  Initial revision
//
//  Revision 1.13  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//

// ***************************************************************************

#endif // BOOST_UNIT_TEST_MONITOR_HPP
