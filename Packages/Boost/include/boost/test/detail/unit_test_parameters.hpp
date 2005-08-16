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
//  Description : storage for unit test framework parameters information
// ***************************************************************************

#ifndef BOOST_TEST_UNIT_TEST_PARAMETERS_HPP_071894GER
#define BOOST_TEST_UNIT_TEST_PARAMETERS_HPP_071894GER

#include <boost/test/detail/global_typedef.hpp>
#include <boost/test/detail/log_level.hpp>

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

namespace boost {

namespace unit_test {

// ************************************************************************** //
// **************                 runtime_config               ************** //
// ************************************************************************** //

namespace runtime_config {

void                    init( int* argc, char** argv );

unit_test::log_level    log_level();
bool                    no_result_code();
unit_test::report_level report_level();
const_string            test_to_run();
bool                    save_pattern();
bool                    show_build_info();
bool                    show_progress();
bool                    catch_sys_errors();
output_format           report_format();
output_format           log_format();
long                    detect_memory_leak();
int                     random_seed();

} // namespace runtime_config

} // namespace unit_test

} // namespace boost

//____________________________________________________________________________//

#include <boost/test/detail/enable_warnings.hpp>

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.3  2005/08/16 11:24:13  spion
//  Import of Boost v. 1.33.0
//
//  Revision 1.21  2005/04/05 06:11:37  rogeeff
//  memory leak allocation point detection\nextra help with _WIN32_WINNT
//
//  Revision 1.20  2005/02/21 10:18:30  rogeeff
//  random cla support
//
//  Revision 1.19  2005/02/20 08:27:06  rogeeff
//  This a major update for Boost.Test framework. See release docs for complete list of fixes/updates
//
// ***************************************************************************

#endif // BOOST_TEST_UNIT_TEST_PARAMETERS_HPP_071894GER
