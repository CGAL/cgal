/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(INTERPRET_PRAGMA_HPP_B1F2315E_C5CE_4ED1_A343_0EF548B7942A_INCLUDED)
#define INTERPRET_PRAGMA_HPP_B1F2315E_C5CE_4ED1_A343_0EF548B7942A_INCLUDED

#include <string>
#include <list>

#include <boost/spirit/core.hpp>
#if SPIRIT_VERSION >= 0x1700
#include <boost/spirit/actor/assign_actor.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#endif // SPIRIT_VERSION >= 0x1700

#include <boost/wave/wave_config.hpp>

#include <boost/wave/util/pattern_parser.hpp>
#include <boost/wave/util/macro_helpers.hpp>

#include <boost/wave/token_ids.hpp>
#include <boost/wave/cpp_exceptions.hpp>
#include <boost/wave/cpp_iteration_context.hpp>
#include <boost/wave/language_support.hpp>

#if !defined(spirit_append_actor)
#if SPIRIT_VERSION >= 0x1700
#define spirit_append_actor(actor) boost::spirit::push_back_a(actor)
#define spirit_assign_actor(actor) boost::spirit::assign_a(actor)
#else
#define spirit_append_actor(actor) boost::spirit::append(actor)
#define spirit_assign_actor(actor) boost::spirit::assign(actor)
#endif // SPIRIT_VERSION >= 0x1700
#endif // !defined(spirit_append_actor)

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace util {

///////////////////////////////////////////////////////////////////////////////
//
//  The function interpret_pragma interprets the given token sequence as the
//  body of a #pragma directive (or parameter to the _Pragma operator) and 
//  executes the actions associated with recognized Wave specific options.
//
///////////////////////////////////////////////////////////////////////////////
template <typename ContextT, typename IteratorT, typename ContainerT>
inline bool 
interpret_pragma(ContextT &ctx, typename ContextT::token_type const &act_token,
    IteratorT it, IteratorT const &end, ContainerT &pending)
{
    typedef typename ContextT::token_type token_type;
    typedef typename token_type::string_type string_type;
    
    using namespace cpplexer;
    if (T_IDENTIFIER == token_id(*it) && "wave" == (*it).get_value()) {
    //  this is a wave specific option, it should have the form:
    //      #pragma wave option(value)
    //  where '(value)' is required only for some pragma directives
    //  all of the given #pragma operators are forwarded to the supplied 
    //  context_policy    
        using namespace boost::spirit;
        token_type option;
        ContainerT values;
        
        if (!parse (++it, end, 
                        (   ch_p(T_IDENTIFIER)
                            [
                                spirit_assign_actor(option)
                            ] 
                        |   pattern_p(KeywordTokenType, TokenTypeMask)
                            [
                                spirit_assign_actor(option)
                            ] 
                        |   pattern_p(OperatorTokenType|AltExtTokenType, 
                                ExtTokenTypeMask)   // and, bit_and etc.
                            [
                                spirit_assign_actor(option)
                            ] 
                        )
                    >> !(   ch_p(T_LEFTPAREN) 
                        >>  lexeme_d[
                                *(anychar_p[spirit_append_actor(values)] - ch_p(T_RIGHTPAREN))
                            ]
                        >>  ch_p(T_RIGHTPAREN)
                        ),
                pattern_p(WhiteSpaceTokenType, TokenTypeMask)).hit)
        {
            return false;
        }
    
    // remove the falsely matched closing parenthesis
        if (values.size() > 0) {
            if (T_RIGHTPAREN == values.back()) {
            typename ContainerT::reverse_iterator rit = values.rbegin();
            
                values.erase((++rit).base());
            }
            else {
                BOOST_WAVE_THROW(preprocess_exception, ill_formed_pragma_option,
                    "missing matching ')'", act_token.get_position());
            }
        }
        
    // decode the option (call the context_policy hook)
        if (!ctx.interpret_pragma(pending, option, values, act_token)) 
        {
        // unknown #pragma option 
        string_type option_str (option.get_value());

            if (values.size() > 0) {
                option_str += "(";
                option_str += impl::as_string(values);
                option_str += ")";
            }
            BOOST_WAVE_THROW(preprocess_exception, ill_formed_pragma_option,
                option_str.c_str(), act_token.get_position());
        }
        return true;
    }
#if BOOST_WAVE_SUPPORT_PRAGMA_ONCE != 0
    else if (T_IDENTIFIER == token_id(*it) && "once" == (*it).get_value()) {
    // #pragma once
        return ctx.add_pragma_once_header(ctx.get_current_filename());
    }
#endif 

    return false;
}

///////////////////////////////////////////////////////////////////////////////
}   // namespace util
}   // namespace wave
}   // namespace boost

#endif // !defined(INTERPRET_PRAGMA_HPP_B1F2315E_C5CE_4ED1_A343_0EF548B7942A_INCLUDED)
