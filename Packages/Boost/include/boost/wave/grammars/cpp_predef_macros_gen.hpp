/*=============================================================================
    A Standard compliant C++ preprocessor

    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(CPP_PREDEF_MACROS_GEN_HPP_CADB6D2C_76A4_4988_83E1_EFFC6902B9A2_INCLUDED)
#define CPP_PREDEF_MACROS_GEN_HPP_CADB6D2C_76A4_4988_83E1_EFFC6902B9A2_INCLUDED

#include <boost/spirit/tree/parse_tree.hpp>

#include <boost/wave/wave_config.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace grammars {

///////////////////////////////////////////////////////////////////////////////
//  
//  store parser_id's of all rules of the predefined_macros_grammar here 
//  for later access
//
///////////////////////////////////////////////////////////////////////////////
struct predefined_macros_grammar_rule_ids {
    std::size_t plain_define_id;       // #define
    std::size_t macro_parameters_id;
    std::size_t macro_definition_id;
};

///////////////////////////////////////////////////////////////////////////////
//  
//  predefined_macros_grammar_gen template class
//
//      This template helps separating the compilation of the 
//      predefined_macros_grammar class from the compilation of the 
//      main pp_iterator. This is done to safe compilation time.
//
//      This class helps parsing command line given macro definitions in a
//      similar way, as macros are parsed by the cpp_grammar class.
//
///////////////////////////////////////////////////////////////////////////////

template <typename LexIteratorT>
struct predefined_macros_grammar_gen
{
    typedef LexIteratorT iterator_type;

//  the parser_id's of all rules of the cpp_grammar are stored here
//  note: these are valid only after the first call to parse_cpp_grammar
    static predefined_macros_grammar_rule_ids rule_ids;

//  parse the cpp_grammar and return the resulting parse tree    
    static boost::spirit::tree_parse_info<iterator_type> 
    parse_predefined_macro (iterator_type const &first, iterator_type const &last);
};

///////////////////////////////////////////////////////////////////////////////
//  definitions of the static members
template <typename LexIteratorT>
predefined_macros_grammar_rule_ids 
    predefined_macros_grammar_gen<LexIteratorT>::rule_ids;

///////////////////////////////////////////////////////////////////////////////
}   // namespace grammars
}   // namespace wave
}   // namespace boost

#endif // !defined(CPP_PREDEF_MACROS_GEN_HPP_CADB6D2C_76A4_4988_83E1_EFFC6902B9A2_INCLUDED)
