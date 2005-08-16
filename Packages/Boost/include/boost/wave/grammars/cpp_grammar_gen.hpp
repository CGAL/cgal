/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(CPP_GRAMMAR_GEN_HPP_80CB8A59_5411_4E45_B406_62531A12FB99_INCLUDED)
#define CPP_GRAMMAR_GEN_HPP_80CB8A59_5411_4E45_B406_62531A12FB99_INCLUDED

#include <boost/spirit/tree/parse_tree.hpp>

#include <boost/wave/wave_config.hpp>
#include <boost/wave/language_support.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace grammars {

///////////////////////////////////////////////////////////////////////////////
//  
//  store parser_id's of all rules of the cpp_grammar here for later access
//
///////////////////////////////////////////////////////////////////////////////
struct cpp_grammar_rule_ids {
    std::size_t pp_statement_id;
    std::size_t include_file_id;       // #include "..."
    std::size_t sysinclude_file_id;    // #include <...>
    std::size_t macroinclude_file_id;  // #include ...
    std::size_t plain_define_id;       // #define
    std::size_t macro_parameters_id;
    std::size_t macro_definition_id;
    std::size_t undefine_id;           // #undef
    std::size_t ifdef_id;              // #ifdef
    std::size_t ifndef_id;             // #ifndef
    std::size_t if_id;                 // #if
    std::size_t elif_id;               // #elif
    std::size_t else_id;               // #else
    std::size_t endif_id;              // #endif
    std::size_t line_id;               // #line
    std::size_t error_id;              // #error
    std::size_t warning_id;            // #warning
    std::size_t pragma_id;             // #pragma
    std::size_t illformed_id;
    std::size_t ppspace_id;
    std::size_t ppqualifiedname_id;
#if BOOST_WAVE_SUPPORT_MS_EXTENSIONS != 0
    std::size_t region_id;             // #region
    std::size_t endregion_id;          // #endregion
#endif
};

///////////////////////////////////////////////////////////////////////////////
//  
//  cpp_grammar_gen template class
//
//      This template helps separating the compilation of the cpp_grammar 
//      class from the compilation of the main pp_iterator. This is done to
//      safe compilation time.
//
///////////////////////////////////////////////////////////////////////////////

template <typename LexIteratorT>
struct cpp_grammar_gen
{
    typedef LexIteratorT                          iterator_type;
    typedef typename LexIteratorT::token_type     token_type;
    typedef typename token_type::position_type    position_type;
    
//  the parser_id's of all rules of the cpp_grammar are stored here
//  note: these are valid only after the first call to parse_cpp_grammar
    static cpp_grammar_rule_ids rule_ids;

//  the actual position of the last matched T_NEWLINE is stored here into the
//  member 'pos_of_newline'
    static position_type pos_of_newline;

//  the found_eof flag is set to true during the parsing, if the directive 
//  under inspection terminates with a T__EOF token
    static bool found_eof;

//  the found_directive contains the token_id of the recognized pp directive
    static boost::wave::token_id found_directive;
        
//  parse the cpp_grammar and return the resulting parse tree    
    static boost::spirit::tree_parse_info<iterator_type> 
    parse_cpp_grammar (iterator_type const &first, iterator_type const &last,
        bool &found_eof_, position_type const &act_pos);
};

///////////////////////////////////////////////////////////////////////////////
//  definitions of the static members
template <typename LexIteratorT>
cpp_grammar_rule_ids 
    cpp_grammar_gen<LexIteratorT>::rule_ids;

template <typename LexIteratorT>
typename LexIteratorT::token_type::position_type 
    cpp_grammar_gen<LexIteratorT>::pos_of_newline;

template <typename LexIteratorT>
bool cpp_grammar_gen<LexIteratorT>::found_eof = false;

template <typename LexIteratorT>
boost::wave::token_id cpp_grammar_gen<LexIteratorT>::found_directive = 
    boost::wave::T_EOF;

///////////////////////////////////////////////////////////////////////////////
}   // namespace grammars
}   // namespace wave
}   // namespace boost

#endif // !defined(CPP_GRAMMAR_GEN_HPP_80CB8A59_5411_4E45_B406_62531A12FB99_INCLUDED)
