/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(CPP_LITERAL_GRAMMAR_GEN_HPP_67794A6C_468A_4AAB_A757_DEDDB182F5A0_INCLUDED)
#define CPP_LITERAL_GRAMMAR_GEN_HPP_67794A6C_468A_4AAB_A757_DEDDB182F5A0_INCLUDED

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace grammars {

///////////////////////////////////////////////////////////////////////////////
//  
//  cpp_intlit_grammar_gen template class
//
//      This template helps separating the compilation of the intlit_grammar 
//      class from the compilation of the expression_grammar. This is done 
//      to safe compilation time.
//
///////////////////////////////////////////////////////////////////////////////
template <typename TokenT>
struct intlit_grammar_gen {

    static unsigned long evaluate(TokenT const &tok, bool &is_unsigned);
};

///////////////////////////////////////////////////////////////////////////////
//  
//  cpp_chlit_grammar_gen template class
//
//      This template helps separating the compilation of the chlit_grammar 
//      class from the compilation of the expression_grammar. This is done 
//      to safe compilation time.
//
///////////////////////////////////////////////////////////////////////////////
template <typename TokenT>
struct chlit_grammar_gen {

    static unsigned int evaluate(TokenT const &tok);
};

///////////////////////////////////////////////////////////////////////////////
}   //  namespace grammars
}   //  namespace wave 
}   //  namespace boost

#endif // !defined(CPP_LITERAL_GRAMMAR_GEN_HPP_67794A6C_468A_4AAB_A757_DEDDB182F5A0_INCLUDED)
