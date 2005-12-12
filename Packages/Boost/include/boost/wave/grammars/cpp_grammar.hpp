/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(CPP_GRAMMAR_HPP_FEAEBC2E_2734_428B_A7CA_85E5A415E23E_INCLUDED)
#define CPP_GRAMMAR_HPP_FEAEBC2E_2734_428B_A7CA_85E5A415E23E_INCLUDED

#include <boost/spirit/core.hpp>
#include <boost/spirit/tree/parse_tree.hpp>
#include <boost/spirit/tree/parse_tree_utils.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <boost/spirit/utility/lists.hpp>

#include <boost/wave/wave_config.hpp>

#if BOOST_WAVE_DUMP_PARSE_TREE != 0
#include <map>
#include <boost/spirit/tree/tree_to_xml.hpp>
#endif

#include <boost/wave/token_ids.hpp>
#include <boost/wave/grammars/cpp_grammar_gen.hpp>
#include <boost/wave/util/pattern_parser.hpp>

#include <boost/wave/cpp_exceptions.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave { 
namespace grammars {

namespace impl {

///////////////////////////////////////////////////////////////////////////////
//
//  store_position 
//
//      The store_position functor extracts the actual file position from the 
//      supplied token.
//
///////////////////////////////////////////////////////////////////////////////

    template <typename PositionT>
    struct store_position {

        store_position(PositionT &pos_) : pos(pos_) {}
        
        template <typename TokenT>
        void operator()(TokenT const &token) const
        {
            pos = token.get_position();
        }
        
        PositionT &pos;
    };

///////////////////////////////////////////////////////////////////////////////
//
//  store_found_eof
//
//      The store_found_eof functor sets a given flag if the T_EOF token was 
//      found during the parsing process
//
///////////////////////////////////////////////////////////////////////////////

    struct store_found_eof {

        store_found_eof(bool &found_eof_) : found_eof(found_eof_) {}
        
        template <typename TokenT>
        void operator()(TokenT const &token) const
        {
            found_eof = true;
        }
        
        bool &found_eof;
    };

///////////////////////////////////////////////////////////////////////////////
//
//  store_found_directive
//
//      The store_found_directive functor stores the token_id of the recognized
//      pp directive
//
///////////////////////////////////////////////////////////////////////////////

    struct store_found_directive {

        store_found_directive(boost::wave::token_id &found_directive_) 
        :   found_directive(found_directive_) {}
        
        template <typename TokenT>
        void operator()(TokenT const &token) const
        {
            found_directive = boost::wave::token_id(token);
        }
        
        boost::wave::token_id &found_directive;
    };

///////////////////////////////////////////////////////////////////////////////
//
//  flush_underlying_parser
//
//      The flush_underlying_parser flushes the underlying
//      multi_pass_iterator during the normal parsing process. This is
//      used at certain points during the parsing process, when it is
//      clear, that no backtracking is needed anymore and the input
//      gathered so far may be discarded.
//
///////////////////////////////////////////////////////////////////////////////
    struct flush_underlying_parser
    :   public boost::spirit::parser<flush_underlying_parser>
    {
        typedef flush_underlying_parser this_t;

        template <typename ScannerT>
        typename boost::spirit::parser_result<this_t, ScannerT>::type
        parse(ScannerT const& scan) const
        {
            scan.first.clear_queue();
            return scan.empty_match();  
        }
    };

    flush_underlying_parser const 
        flush_underlying_parser_p = flush_underlying_parser();

}   // anonymous namespace

///////////////////////////////////////////////////////////////////////////////
//  define, whether the rule's should generate some debug output
#define TRACE_CPP_GRAMMAR \
    bool(BOOST_SPIRIT_DEBUG_FLAGS_CPP & BOOST_SPIRIT_DEBUG_FLAGS_CPP_GRAMMAR) \
    /**/

///////////////////////////////////////////////////////////////////////////////
// Encapsulation of the C++ preprocessor grammar.
template <typename PositionT>
struct cpp_grammar : 
    public boost::spirit::grammar<cpp_grammar<PositionT> >
{
    typedef cpp_grammar<PositionT>          grammar_t;
    typedef impl::store_position<PositionT> store_pos_t;
    typedef impl::store_found_eof           store_found_eof_t;
    typedef impl::store_found_directive     store_found_directive_t;
    
    template <typename ScannerT>
    struct definition
    {
    // non-parse_tree generating rule type
        typedef typename ScannerT::iteration_policy_t iteration_policy_t;
        typedef boost::spirit::match_policy match_policy_t;
        typedef typename ScannerT::action_policy_t action_policy_t;
        typedef 
            boost::spirit::scanner_policies<
                iteration_policy_t, match_policy_t, action_policy_t> 
            policies_t;
        typedef 
            boost::spirit::scanner<typename ScannerT::iterator_t, policies_t> 
            non_tree_scanner_t;
        typedef boost::spirit::rule<non_tree_scanner_t> no_tree_rule_t;

    // 'normal' (parse_tree generating) rule type
        typedef boost::spirit::rule<ScannerT> rule_t;

        rule_t pp_statement;
        rule_t include_file, system_include_file, macro_include_file;
        rule_t plain_define, macro_definition, macro_parameters;
        rule_t undefine;
        rule_t ppifdef, ppifndef, ppif, ppelse, ppelif, ppendif;
        rule_t ppline; 
        rule_t pperror;
        rule_t ppwarning;
        rule_t pppragma;
        rule_t illformed;
        rule_t ppqualifiedname;
        rule_t eol_tokens;
        no_tree_rule_t ppsp;
#if BOOST_WAVE_SUPPORT_MS_EXTENSIONS != 0
        rule_t ppregion;
        rule_t ppendregion;
#endif

        definition(cpp_grammar const &self) 
        {
        // import the spirit and cpplexer namespaces here
            using namespace boost::spirit;
            using namespace boost::wave;
            using namespace boost::wave::util;

        // save the rule id's for later use
            self.rule_ids.pp_statement_id = pp_statement.id().to_long();
            self.rule_ids.include_file_id = include_file.id().to_long();
            self.rule_ids.sysinclude_file_id = system_include_file.id().to_long();
            self.rule_ids.macroinclude_file_id = macro_include_file.id().to_long();
            self.rule_ids.plain_define_id = plain_define.id().to_long();
            self.rule_ids.macro_parameters_id = macro_parameters.id().to_long();
            self.rule_ids.macro_definition_id = macro_definition.id().to_long();
            self.rule_ids.undefine_id = undefine.id().to_long();
            self.rule_ids.ifdef_id = ppifdef.id().to_long();
            self.rule_ids.ifndef_id = ppifndef.id().to_long();
            self.rule_ids.if_id = ppif.id().to_long();
            self.rule_ids.elif_id = ppelif.id().to_long();
            self.rule_ids.else_id = ppelse.id().to_long();
            self.rule_ids.endif_id = ppendif.id().to_long();
            self.rule_ids.line_id = ppline.id().to_long();
            self.rule_ids.error_id = pperror.id().to_long();
            self.rule_ids.warning_id = ppwarning.id().to_long();
            self.rule_ids.pragma_id = pppragma.id().to_long();
            self.rule_ids.illformed_id = illformed.id().to_long();
            self.rule_ids.ppspace_id = ppsp.id().to_long();
            self.rule_ids.ppqualifiedname_id = ppqualifiedname.id().to_long();
#if BOOST_WAVE_SUPPORT_MS_EXTENSIONS != 0
            self.rule_ids.region_id = ppregion.id().to_long();
            self.rule_ids.endregion_id = ppendregion.id().to_long();
#endif

#if BOOST_WAVE_DUMP_PARSE_TREE != 0
            self.map_rule_id_to_name.init_rule_id_to_name_map(self);
#endif 

        // recognizes preprocessor directives only

        // C++ standard 16.1: A preprocessing directive consists of a sequence 
        // of preprocessing tokens. The first token in the sequence is # 
        // preprocessing token that is either the first character in the source 
        // file (optionally after white space containing no new-line 
        // characters) or that follows white space containing at least one 
        // new-line character. The last token in the sequence is the first 
        // new-line character that follows the first token in the sequence.

            pp_statement
                =   (   include_file
                    |   system_include_file
                    |   macro_include_file
                    |   plain_define
                    |   undefine
                    |   ppifdef
                    |   ppifndef
                    |   ppif
                    |   ppelse
                    |   ppelif
                    |   ppendif
                    |   ppline
                    |   pperror
                    |   ppwarning
                    |   pppragma
                    |   illformed
#if BOOST_WAVE_SUPPORT_MS_EXTENSIONS != 0
                    |   ppregion
                    |   ppendregion
#endif
                    )
                    >> eol_tokens
//  In parser debug mode it is useful not to flush the underlying stream
//  to allow its investigation in the debugger and to see the correct
//  output in the printed debug log..
//  Note: this may break the parser, though.
#if !(defined(BOOST_SPIRIT_DEBUG) && \
      (BOOST_SPIRIT_DEBUG_FLAGS_CPP & BOOST_SPIRIT_DEBUG_FLAGS_CPP_GRAMMAR) \
     )
                    >>  impl::flush_underlying_parser_p
#endif // !(defined(BOOST_SPIRIT_DEBUG) &&
                ;

        // #include ...
            include_file            // include "..."
                =   ch_p(T_PP_QHEADER) 
                    [ store_found_directive_t(self.found_directive) ]
#if BOOST_WAVE_SUPPORT_INCLUDE_NEXT != 0
                |   ch_p(T_PP_QHEADER_NEXT)
                    [ store_found_directive_t(self.found_directive) ]
#endif 
                ;

            system_include_file     // include <...>
                =   ch_p(T_PP_HHEADER) 
                    [ store_found_directive_t(self.found_directive) ]
#if BOOST_WAVE_SUPPORT_INCLUDE_NEXT != 0
                |   ch_p(T_PP_HHEADER_NEXT)
                    [ store_found_directive_t(self.found_directive) ]
#endif 
                ;

            macro_include_file      // include ...anything else...
                =   no_node_d
                    [
                        ch_p(T_PP_INCLUDE)
                        [ store_found_directive_t(self.found_directive) ]
#if BOOST_WAVE_SUPPORT_INCLUDE_NEXT != 0
                    |   ch_p(T_PP_INCLUDE_NEXT)
                        [ store_found_directive_t(self.found_directive) ]
#endif
                    ]
                    >> *(   anychar_p -
                            (ch_p(T_NEWLINE) | ch_p(T_CPPCOMMENT) | ch_p(T_EOF)) 
                        )
                ;

        // #define FOO foo (with optional parameters)
            plain_define
                =   no_node_d
                    [
                        ch_p(T_PP_DEFINE) 
                        [ store_found_directive_t(self.found_directive) ]
                        >> +ppsp
                    ]
                    >>  (   ch_p(T_IDENTIFIER) 
                        |   pattern_p(KeywordTokenType, TokenTypeMask)
                        |   pattern_p(OperatorTokenType|AltExtTokenType, 
                                ExtTokenTypeMask)   // and, bit_and etc.
                        )
                    >>  (   (   no_node_d[eps_p(ch_p(T_LEFTPAREN))]
                                >>  macro_parameters
                                >> !macro_definition
                            )
                        |  !(   no_node_d[+ppsp] 
                                >>  macro_definition
                            )
                        )
                ;

        // parameter list
        // normal C++ mode
            macro_parameters
                =   confix_p(
                        no_node_d[ch_p(T_LEFTPAREN) >> *ppsp],
                        !list_p(
                            (   ch_p(T_IDENTIFIER) 
                            |   pattern_p(KeywordTokenType, TokenTypeMask)
                            |   pattern_p(OperatorTokenType|AltExtTokenType, 
                                    ExtTokenTypeMask)   // and, bit_and etc.
#if BOOST_WAVE_SUPPORT_VARIADICS_PLACEMARKERS != 0
                            |   ch_p(T_ELLIPSIS)
#endif
                            ), 
                            no_node_d[*ppsp >> ch_p(T_COMMA) >> *ppsp]
                        ),
                        no_node_d[*ppsp >> ch_p(T_RIGHTPAREN)]
                    )
                ;
            
        // macro body (anything left until eol)
            macro_definition
                =   no_node_d[*ppsp]
                    >> *(   anychar_p -
                            (ch_p(T_NEWLINE) | ch_p(T_CPPCOMMENT) | ch_p(T_EOF)) 
                        )
                ;

        // #undef FOO 
            undefine
                =   no_node_d
                    [
                        ch_p(T_PP_UNDEF) 
                        [ store_found_directive_t(self.found_directive) ]
                        >> +ppsp
                    ]
                    >>  (   ch_p(T_IDENTIFIER) 
                        |   pattern_p(KeywordTokenType, TokenTypeMask)
                        |   pattern_p(OperatorTokenType|AltExtTokenType, 
                                ExtTokenTypeMask)   // and, bit_and etc.
                        )
                ;

        // #ifdef et.al.
            ppifdef
                =   no_node_d
                    [
                        ch_p(T_PP_IFDEF) 
                        [ store_found_directive_t(self.found_directive) ]
                        >> +ppsp
                    ]
                    >>  ppqualifiedname
                ;

            ppifndef
                =   no_node_d
                    [
                        ch_p(T_PP_IFNDEF) 
                        [ store_found_directive_t(self.found_directive) ]
                        >> +ppsp
                    ]
                    >>  ppqualifiedname
                ;

            ppif
                =   no_node_d
                    [
                        ch_p(T_PP_IF) 
                        [ store_found_directive_t(self.found_directive) ]
                        >> *ppsp
                    ]
                    >> +(   anychar_p -
                            (ch_p(T_NEWLINE) | ch_p(T_CPPCOMMENT) | ch_p(T_EOF)) 
                        )
                ;

            ppelse
                =   no_node_d
                    [
                        ch_p(T_PP_ELSE)
                        [ store_found_directive_t(self.found_directive) ]
                    ]
                ;

            ppelif
                =   no_node_d
                    [
                        ch_p(T_PP_ELIF) 
                        [ store_found_directive_t(self.found_directive) ]
                        >> *ppsp
                    ]
                    >> +(   anychar_p -
                            (ch_p(T_NEWLINE) | ch_p(T_CPPCOMMENT) | ch_p(T_EOF)) 
                        )
                ;

            ppendif
                =   no_node_d
                    [
                        ch_p(T_PP_ENDIF)
                        [ store_found_directive_t(self.found_directive) ]
                    ]
                ;

        // #line ...
            ppline 
                =   no_node_d
                    [
                        ch_p(T_PP_LINE) 
                        [ store_found_directive_t(self.found_directive) ]
                        >> *ppsp
                    ]
                    >> +(   anychar_p -
                            (ch_p(T_NEWLINE) | ch_p(T_CPPCOMMENT) | ch_p(T_EOF)) 
                        )
                ;
                
#if BOOST_WAVE_SUPPORT_MS_EXTENSIONS != 0
        // #region ...
            ppregion
                =   no_node_d
                    [
                        ch_p(T_MSEXT_PP_REGION) 
                        [ store_found_directive_t(self.found_directive) ]
                        >> +ppsp
                    ]
                    >>  ppqualifiedname
                ;

        // #endregion
            ppendregion
                =   no_node_d
                    [
                        ch_p(T_MSEXT_PP_ENDREGION) 
                        [ store_found_directive_t(self.found_directive) ]
                    ]
                ;
#endif

        // # something else (ill formed preprocessor directive)
            illformed           // for error reporting
                =   no_node_d
                    [
                        pattern_p(T_POUND, MainTokenMask) 
                        >> *ppsp
                    ]
                    >>  (   anychar_p -
                            (ch_p(T_NEWLINE) | ch_p(T_CPPCOMMENT) | ch_p(T_EOF)) 
                        )
                    >>  no_node_d
                        [
                           *(   anychar_p -
                                (ch_p(T_NEWLINE) | ch_p(T_CPPCOMMENT) | ch_p(T_EOF)) 
                            )
                        ]
                ;

        // #error
            pperror
                =   no_node_d
                    [
                        ch_p(T_PP_ERROR) 
                        [ store_found_directive_t(self.found_directive) ]
                        >> *ppsp
                    ]
                    >> *(   anychar_p -
                            (ch_p(T_NEWLINE) | ch_p(T_CPPCOMMENT) | ch_p(T_EOF)) 
                        )
                ;

        // #warning
            ppwarning
                =   no_node_d
                    [
                        ch_p(T_PP_WARNING) 
                        [ store_found_directive_t(self.found_directive) ]
                        >> *ppsp
                    ]
                    >> *(   anychar_p -
                            (ch_p(T_NEWLINE) | ch_p(T_CPPCOMMENT) | ch_p(T_EOF)) 
                        )
                ;

        // #pragma ...
            pppragma
                =   no_node_d
                    [
                        ch_p(T_PP_PRAGMA)
                        [ store_found_directive_t(self.found_directive) ]
                    ] 
                    >> *(   anychar_p -
                            (ch_p(T_NEWLINE) | ch_p(T_CPPCOMMENT) | ch_p(T_EOF)) 
                        )
                ;

            ppqualifiedname
                =   no_node_d[*ppsp]
                    >>  (   ch_p(T_IDENTIFIER) 
                        |   pattern_p(KeywordTokenType, TokenTypeMask)
                        |   pattern_p(OperatorTokenType|AltExtTokenType, 
                                ExtTokenTypeMask)   // and, bit_and etc.
                        ) 
                ;

        // auxiliary helper rules
            ppsp     // valid space in a line with a preprocessor directive
                =   ch_p(T_SPACE) | ch_p(T_CCOMMENT)
                ;

        // end of line tokens
            eol_tokens 
                =   no_node_d
                    [
                        *ppsp 
                        >>  (   ch_p(T_NEWLINE)
                                [ store_pos_t(self.pos_of_newline) ]
                            |   ch_p(T_CPPCOMMENT)
                                [ store_pos_t(self.pos_of_newline) ]
                            |   ch_p(T_EOF)
                                [ store_pos_t(self.pos_of_newline) ]
                                [ store_found_eof_t(self.found_eof) ]
                            )
                    ]
                ;

            BOOST_SPIRIT_DEBUG_TRACE_RULE(pp_statement, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(include_file, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(system_include_file, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(macro_include_file, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(plain_define, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(macro_definition, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(macro_parameters, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(undefine, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppifdef, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppifndef, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppif, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppelse, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppelif, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppendif, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppline, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(pperror, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppwarning, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(illformed, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppsp, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppqualifiedname, TRACE_CPP_GRAMMAR);
#if BOOST_WAVE_SUPPORT_MS_EXTENSIONS != 0
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppregion, TRACE_CPP_GRAMMAR);
            BOOST_SPIRIT_DEBUG_TRACE_RULE(ppendregion, TRACE_CPP_GRAMMAR);
#endif
        }

    // start rule of this grammar
        rule_t const& start() const
        { return pp_statement; }
    };

    cpp_grammar_rule_ids &rule_ids;
    PositionT &pos_of_newline;
    bool &found_eof;
    boost::wave::token_id &found_directive;
    
    cpp_grammar(cpp_grammar_rule_ids &rule_ids_, PositionT &pos_of_newline_,
            bool &found_eof_, boost::wave::token_id &found_directive_)
    :   rule_ids(rule_ids_), pos_of_newline(pos_of_newline_), 
        found_eof(found_eof_), found_directive(found_directive_)
    { 
        BOOST_SPIRIT_DEBUG_TRACE_GRAMMAR_NAME(*this, "cpp_grammar", 
            TRACE_CPP_GRAMMAR); 
    }

#if BOOST_WAVE_DUMP_PARSE_TREE != 0
// helper function and data to get readable names of the rules known to us
    struct map_ruleid_to_name :
        public std::map<boost::spirit::parser_id, std::string> 
    {
        typedef std::map<boost::spirit::parser_id, std::string> base_t;

        void init_rule_id_to_name_map(cpp_grammar const &self)
        {
            struct {
                int parser_id;
                char const *rule_name;
            } 
            init_ruleid_name_map[] = {
                { self.rule_ids.pp_statement_id, "pp_statement" },
                { self.rule_ids.include_file_id, "include_file" },
                { self.rule_ids.sysinclude_file_id, "system_include_file" },
                { self.rule_ids.macroinclude_file_id, "macro_include_file" },
                { self.rule_ids.plain_define_id, "plain_define" },
                { self.rule_ids.macro_parameters_id, "macro_parameters" },
                { self.rule_ids.macro_definition_id, "macro_definition" },
                { self.rule_ids.undefine_id, "undefine" },
                { self.rule_ids.ifdef_id, "ppifdef" },
                { self.rule_ids.ifndef_id, "ppifndef" },
                { self.rule_ids.if_id, "ppif" },
                { self.rule_ids.elif_id, "ppelif" },
                { self.rule_ids.else_id, "ppelse" },
                { self.rule_ids.endif_id, "ppendif" },
                { self.rule_ids.line_id, "ppline" },
                { self.rule_ids.error_id, "pperror" },
                { self.rule_ids.warning_id, "ppwarning" },
                { self.rule_ids.pragma_id, "pppragma" },
                { self.rule_ids.illformed_id, "illformed" },
                { self.rule_ids.ppqualifiedname_id, "ppqualifiedname" },
#if BOOST_WAVE_SUPPORT_MS_EXTENSIONS != 0
                { self.rule_ids.region_id, "ppregion" },
                { self.rule_ids.endregion_id, "ppendregion" },
#endif
                { 0 }
            };

        // initialize parser_id to rule_name map
            for (int i = 0; 0 != init_ruleid_name_map[i].parser_id; ++i)
                base_t::insert(base_t::value_type(
                    boost::spirit::parser_id(init_ruleid_name_map[i].parser_id), 
                    std::string(init_ruleid_name_map[i].rule_name))
                );
        }
    };
    mutable map_ruleid_to_name map_rule_id_to_name;
#endif // WAVE_DUMP_PARSE_TREE != 0
};

///////////////////////////////////////////////////////////////////////////////
#undef TRACE_CPP_GRAMMAR

///////////////////////////////////////////////////////////////////////////////
//  
//  The following parse function is defined here, to allow the separation of 
//  the compilation of the cpp_grammar from the function using it.
//  
///////////////////////////////////////////////////////////////////////////////

#if BOOST_WAVE_SEPARATE_GRAMMAR_INSTANTIATION != 0
#define BOOST_WAVE_GRAMMAR_GEN_INLINE
#else
#define BOOST_WAVE_GRAMMAR_GEN_INLINE inline
#endif 

namespace {

    char const *get_directivename(boost::wave::token_id id)
    {
        using namespace boost::wave;
        switch (static_cast<unsigned int>(id)) {
        case T_PP_QHEADER:
        case T_PP_HHEADER:
        case T_PP_INCLUDE:  return "#include";
        case T_PP_DEFINE:   return "#define";
        case T_PP_UNDEF:    return "#undef";
        case T_PP_IFDEF:    return "#ifdef";
        case T_PP_IFNDEF:   return "#ifndef";
        case T_PP_IF:       return "#if";
        case T_PP_ELSE:     return "#else";
        case T_PP_ELIF:     return "#elif";
        case T_PP_ENDIF:    return "#endif";
        case T_PP_LINE:     return "#line";
        case T_PP_ERROR:    return "#error";
        case T_PP_WARNING:  return "#warning";
        case T_PP_PRAGMA:   return "#pragma";
        default:
            return "#unknown directive";
        }
    }
}

template <typename LexIteratorT>
BOOST_WAVE_GRAMMAR_GEN_INLINE 
boost::spirit::tree_parse_info<LexIteratorT>
cpp_grammar_gen<LexIteratorT>::parse_cpp_grammar (
    LexIteratorT const &first, LexIteratorT const &last,
    bool &found_eof_, position_type const &act_pos)
{
    using namespace boost::spirit;
    using namespace boost::wave;
    
    pos_of_newline = position_type();  // reset position
    found_eof = false;              // reset flag
    found_directive = T_EOF;        // reset found directive
    
    cpp_grammar<position_type> g(
        rule_ids, pos_of_newline, found_eof, found_directive);
    
    tree_parse_info<LexIteratorT> hit = pt_parse (first, last, g);
    
#if BOOST_WAVE_DUMP_PARSE_TREE != 0
    if (hit.match) {
        tree_to_xml (BOOST_WAVE_DUMP_PARSE_TREE_OUT, hit.trees, "", 
            g.map_rule_id_to_name, &token_type::get_token_id, 
            &token_type::get_token_value);
    }
#endif

    if (!hit.match && found_directive != T_EOF) {
    // recognized invalid directive
    std::string directive = get_directivename(found_directive);
    
        BOOST_WAVE_THROW(preprocess_exception, ill_formed_directive, 
            directive.c_str(), act_pos);
    }

    found_eof_ = found_eof;
    return hit;
}

#undef BOOST_WAVE_GRAMMAR_GEN_INLINE

///////////////////////////////////////////////////////////////////////////////
}   // namespace grammars
}   // namespace wave
}   // namespace boost 

#endif // !defined(CPP_GRAMMAR_HPP_FEAEBC2E_2734_428B_A7CA_85E5A415E23E_INCLUDED)
