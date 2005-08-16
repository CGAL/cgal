//  (C) Copyright Gennadiy Rozental 2001-2005.
//  (C) Copyright Beman Dawes 1995-2001.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : main function implementation for Program Executon Monitor
// ***************************************************************************

#ifndef BOOST_TEST_CPP_MAIN_IPP_012205GER
#define BOOST_TEST_CPP_MAIN_IPP_012205GER

// Boost.Test
#include <boost/test/execution_monitor.hpp>
#include <boost/test/detail/config.hpp>
#include <boost/test/utils/basic_cstring/io.hpp>

// Boost
#include <boost/cstdlib.hpp>    // for exit codes
#include <boost/config.hpp>     // for workarounds

// STL
#include <iostream>
#include <cstdlib>      // std::getenv

#include <boost/test/detail/suppress_warnings.hpp>

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

#ifdef BOOST_NO_STDC_NAMESPACE
namespace std { using ::getenv; }
#endif

int cpp_main( int argc, char* argv[] );  // prototype for user's cpp_main()

namespace {

struct cpp_main_caller {
    cpp_main_caller( int argc, char** argv ) : m_argc( argc ), m_argv( argv ) {}
    
    int operator()() { return cpp_main( m_argc, m_argv ); }
  
private:
    // Data members    
    int      m_argc;
    char**   m_argv;
};

} // local namespace

// ************************************************************************** //
// **************                   cpp main                   ************** //
// ************************************************************************** //

int BOOST_TEST_CALL_DECL main( int argc, char* argv[] )
{
    int result;

    boost::unit_test::const_string p( std::getenv( "BOOST_TEST_CATCH_SYSTEM_ERRORS" ) );
    bool catch_system_errors = p != "no";
        
    try {
        ::boost::execution_monitor ex_mon;
        result = ex_mon.execute( ::boost::unit_test::callback0<int>( cpp_main_caller( argc, argv ) ), catch_system_errors );
        
        if( result == 0 )
            result = ::boost::exit_success;
        else if( result != ::boost::exit_success ) {
            std::cout << "\n**** error return code: " << result << std::endl;
            result = ::boost::exit_failure;
        }
    }
    catch( ::boost::execution_exception const& exex ) {
        std::cout << "\n**** exception(" << exex.code() << "): " << exex.what() << std::endl;
        result = ::boost::exit_exception_failure;
    }
    
    if( result != ::boost::exit_success ) {
        std::cerr << "******** errors detected; see standard output for details ********" << std::endl;
    }
    else {
        //  Some prefer a confirming message when all is well, while others don't
        //  like the clutter.  Use an environment variable to avoid command
        //  line argument modifications; for use in production programs
        //  that's a no-no in some organizations.
        ::boost::unit_test::const_string p( std::getenv( "BOOST_PRG_MON_CONFIRM" ) );
        if( p != "no" ) { 
            std::cerr << std::flush << "no errors detected" << std::endl; 
        }
    }

    return result;
}

//____________________________________________________________________________//

#include <boost/test/detail/enable_warnings.hpp>

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2005/08/16 11:24:13  spion
//  Initial revision
//
//  Revision 1.5  2005/02/20 08:27:07  rogeeff
//  This a major update for Boost.Test framework. See release docs for complete list of fixes/updates
//
//  Revision 1.4  2005/02/01 06:40:07  rogeeff
//  copyright update
//  old log entries removed
//  minor stilistic changes
//  depricated tools removed
//
//  Revision 1.3  2005/01/31 07:50:06  rogeeff
//  cdecl portability fix
//
//  Revision 1.2  2005/01/31 06:01:50  rogeeff
//  BOOST_TEST_CALL_DECL correctness fixes
//
//  Revision 1.1  2005/01/22 19:22:12  rogeeff
//  implementation moved into headers section to eliminate dependency of included/minimal component on src directory
//
// ***************************************************************************

#endif // BOOST_TEST_CPP_MAIN_IPP_012205GER
