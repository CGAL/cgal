/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    Re2C based C++ lexer
    
    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(CPP_RE_HPP_B76C4F5E_63E9_4B8A_9975_EC32FA6BF027_INCLUDED)
#define CPP_RE_HPP_B76C4F5E_63E9_4B8A_9975_EC32FA6BF027_INCLUDED

#include <boost/wave/token_ids.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace cpplexer {
namespace re2clex {

///////////////////////////////////////////////////////////////////////////////
//  The scanner function to call whenever a new token is requested
boost::wave::token_id scan(Scanner *s);

///////////////////////////////////////////////////////////////////////////////
}   // namespace re2clex
}   // namespace cpplexer
}   // namespace wave
}   // namespace boost

#endif // !defined(CPP_RE_HPP_B76C4F5E_63E9_4B8A_9975_EC32FA6BF027_INCLUDED)
