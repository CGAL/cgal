/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    Definition of the abstract lexer interface
    
    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(CPP_LEX_INTERFACE_HPP_E83F52A4_90AC_4FBE_A9A7_B65F7F94C497_INCLUDED)
#define CPP_LEX_INTERFACE_HPP_E83F52A4_90AC_4FBE_A9A7_B65F7F94C497_INCLUDED

#include <boost/wave/util/file_position.hpp>
#include <boost/wave/language_support.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace cpplexer {

///////////////////////////////////////////////////////////////////////////////
//  
//  new_lexer_gen: generates a new instance of the required C++ lexer
//
///////////////////////////////////////////////////////////////////////////////
template <typename TokenT> struct lex_input_interface; 

template <
    typename IteratorT, 
    typename PositionT = boost::wave::util::file_position_type
>
struct new_lexer_gen
{
//  The NewLexer function allows the opaque generation of a new lexer object.
//  It is coupled to the token type to allow to decouple the lexer/token 
//  configurations at compile time.
    static lex_input_interface<lex_token<PositionT> > *
    new_lexer(IteratorT const &first, IteratorT const &last, 
        PositionT const &pos, boost::wave::language_support language);
};

///////////////////////////////////////////////////////////////////////////////
//
//  The lex_input_interface decouples the lex_iterator_shim from the actual 
//  lexer. This is done to allow compile time reduction.
//  Thanks to JCAB for having this idea.
//
///////////////////////////////////////////////////////////////////////////////

template <typename TokenT>
struct lex_input_interface 
{
    typedef typename TokenT::position_type position_type;
    
    virtual TokenT get() = 0;
    virtual void set_position(position_type const &pos) = 0;

    virtual ~lex_input_interface() {}
    
//  The new_lexer function allows the opaque generation of a new lexer object.
//  It is coupled to the token type to allow to distinguish different 
//  lexer/token configurations at compile time.
    template <typename IteratorT>
    static lex_input_interface *
    new_lexer(IteratorT const &first, IteratorT const &last, 
        position_type const &pos, boost::wave::language_support language)
    { 
        return new_lexer_gen<IteratorT, position_type>::new_lexer (first, last, 
            pos, language); 
    }
};

///////////////////////////////////////////////////////////////////////////////
}   // namespace cpplexer
}   // namespace wave
}   // namespace boost 

#endif // !defined(CPP_LEX_INTERFACE_HPP_E83F52A4_90AC_4FBE_A9A7_B65F7F94C497_INCLUDED)
