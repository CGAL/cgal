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
//  Description : storage for unit test framework parameters information
// ***************************************************************************

#ifndef BOOST_UNIT_TEST_PARAMETERS_HPP
#define BOOST_UNIT_TEST_PARAMETERS_HPP

// STL
#include <string>   // std::string

#include <boost/test/detail/unit_test_config.hpp>

namespace boost {

namespace unit_test_framework {

// framework parameters and there corresponding command-line arguments
c_string_literal const LOG_LEVEL         = "BOOST_TEST_LOG_LEVEL";              // --log_level
c_string_literal const NO_RESULT_CODE    = "BOOST_TEST_RESULT_CODE";            // --result_code
c_string_literal const REPORT_LEVEL      = "BOOST_TEST_REPORT_LEVEL";           // --report_level
c_string_literal const TESTS_TO_RUN      = "BOOST_TESTS_TO_RUN";                // --run_test
c_string_literal const SAVE_TEST_PATTERN = "BOOST_TEST_SAVE_PATTERN";           // --save_pattern
c_string_literal const BUILD_INFO        = "BOOST_TEST_BUILD_INFO";             // --build_info
c_string_literal const CATCH_SYS_ERRORS  = "BOOST_TEST_CATCH_SYSTEM_ERRORS";    // --catch_system_errors
c_string_literal const REPORT_FORMAT     = "BOOST_TEST_REPORT_FORMAT";          // --report_format
c_string_literal const LOG_FORMAT        = "BOOST_TEST_LOG_FORMAT";             // --log_format
c_string_literal const OUTPUT_FORMAT     = "BOOST_TEST_OUTPUT_FORMAT";          // --output_format

enum report_level                             { CONFIRMATION_REPORT, SHORT_REPORT, DETAILED_REPORT, NO_REPORT, UNDEF_REPORT };
c_string_literal const report_level_names[] = { "confirm"          , "short"     , "detailed"     , "no"     };

enum output_format { HRF /* human readable format */, XML /* XML */ };

std::string retrieve_framework_parameter( c_string_literal parameter_name_, int* argc_, char** argv_ );

} // namespace unit_test_framework

} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/05/23 10:51:39  spion
//  Initial revision
//
//  Revision 1.12  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//

// ***************************************************************************

#endif // BOOST_UNIT_TEST_CONFIG_HPP
