//  Copyright (C) Christof Meerwald 2003
//  Copyright (C) Dan Watkins 2003
//
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Digital Mars C++ compiler setup:
#define BOOST_COMPILER __DMC_VERSION_STRING__

#define BOOST_HAS_LONG_LONG
#define BOOST_HAS_PRAGMA_ONCE

#define BOOST_FUNCTION_SCOPE_USING_DECLARATION_BREAKS_ADL
#define BOOST_NO_OPERATORS_IN_NAMESPACE
#define BOOST_NO_SFINAE
#define BOOST_NO_TEMPLATE_TEMPLATES
#define BOOST_NO_USING_TEMPLATE
#define BOOST_NEEDS_TOKEN_PASTING_OP_FOR_TOKENS_JUXTAPOSING
#define BOOST_NO_ARRAY_TYPE_SPECIALIZATIONS

// check for exception handling support:
#ifndef _CPPUNWIND
#  define BOOST_NO_EXCEPTIONS
#endif

#if (__DMC__ < 0x833)
#  if defined(BOOST_ASSERT_CONFIG)
#     error "Unknown compiler version - please run the configure tests and report the results"
#  endif
#endif
