/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    Global application configuration
    
    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(BOOST_SPIRIT_PATTERN_PARSER_HPP)
#define BOOST_SPIRIT_PATTERN_PARSER_HPP

#include <boost/spirit/core/primitives/primitives.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace util {

    ///////////////////////////////////////////////////////////////////////////
    //
    //  pattern_and class
    //
    ///////////////////////////////////////////////////////////////////////////
    template <typename CharT = char>
    struct pattern_and : public boost::spirit::char_parser<pattern_and<CharT> >
    {
        pattern_and(CharT pattern_, unsigned long pattern_mask_ = 0UL)
        :   pattern(pattern_), 
            pattern_mask((0UL != pattern_mask_) ? pattern_mask_ : pattern_)
        {}

        template <typename T>
        bool test(T pattern_) const
        { return (pattern_ & pattern_mask) == pattern; }

        CharT         pattern;
        unsigned long pattern_mask;
    };

    template <typename CharT>
    inline pattern_and<CharT>
    pattern_p(CharT pattern, unsigned long pattern_mask = 0UL)
    { return pattern_and<CharT>(pattern, pattern_mask); }


///////////////////////////////////////////////////////////////////////////////
}   // namespace util
}   // namespace wave
}   // namespace boost

#endif // defined(BOOST_SPIRIT_PATTERN_PARSER_HPP)
