// -- lambda.hpp -- Boost Lambda Library -----------------------------------
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
// For more information, see http://lambda.cs.utu.fi 

#ifndef BOOST_LAMBDA_LAMBDA_HPP
#define BOOST_LAMBDA_LAMBDA_HPP


#include "boost/lambda/core.hpp"

#ifdef BOOST_NO_FDECL_TEMPLATES_AS_TEMPLATE_TEMPLATE_PARAMS
#include <istream>
#include <ostream>
#endif

#include "boost/lambda/detail/operator_actions.hpp"
#include "boost/lambda/detail/operator_lambda_func_base.hpp"
#include "boost/lambda/detail/operator_return_type_traits.hpp"


#include "boost/lambda/detail/operators.hpp"

#ifndef BOOST_LAMBDA_FAILS_IN_TEMPLATE_KEYWORD_AFTER_SCOPE_OPER
// sorry, member ptr does not work with gcc2.95
#include "boost/lambda/detail/member_ptr.hpp"
#endif


#endif
