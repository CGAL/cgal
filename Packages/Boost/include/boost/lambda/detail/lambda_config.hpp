// Boost Lambda Library - lambda_config.hpp ------------------------------

// Copyright (C) 1999, 2000 Jaakko Järvi (jaakko.jarvi@cs.utu.fi)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies. 
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice 
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty, 
// and with no claim as to its suitability for any purpose.
//
// For more information, see www.boost.org

// ---------------------------------------------------------------

#ifndef BOOST_LAMBDA_LAMBDA_CONFIG_HPP
#define BOOST_LAMBDA_LAMBDA_CONFIG_HPP

// add to boost/config.hpp
// for now


# if defined __GNUC__
#   if (__GNUC__ == 2 && __GNUC_MINOR__ <= 97) 
#define BOOST_NO_TEMPLATED_STREAMS
#define BOOST_LAMBDA_INCORRECT_BIND_OVERLOADING
#endif
#   if (__GNUC__ == 2 && __GNUC_MINOR__ <= 95) 
#define BOOST_LAMBDA_FAILS_IN_TEMPLATE_KEYWORD_AFTER_SCOPE_OPER
#endif
#endif  // __GNUC__
 

#if defined __KCC

#define BOOST_NO_FDECL_TEMPLATES_AS_TEMPLATE_TEMPLATE_PARAMS

#endif  // __KCC

#endif







