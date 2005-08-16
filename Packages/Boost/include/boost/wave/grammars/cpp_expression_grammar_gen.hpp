/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(CPP_EXPRESSION_GRAMMAR_GEN_HPP_42399258_6CDC_4101_863D_5C7D95B5A6CA_INCLUDED)
#define CPP_EXPRESSION_GRAMMAR_GEN_HPP_42399258_6CDC_4101_863D_5C7D95B5A6CA_INCLUDED

#include <list>
#include <boost/pool/pool_alloc.hpp>

#include <boost/wave/cpp_iteration_context.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace grammars {

///////////////////////////////////////////////////////////////////////////////
//  
//  expression_grammar_gen template class
//
//      This template helps separating the compilation of the 
//      expression_grammar class from the compilation of the main 
//      pp_iterator. This is done to safe compilation time.
//
///////////////////////////////////////////////////////////////////////////////

template <typename TokenT>
struct expression_grammar_gen {

    typedef TokenT token_type;
    typedef std::list<token_type, boost::fast_pool_allocator<token_type> >
        token_sequence_type;
        
    static bool evaluate(
        typename token_sequence_type::const_iterator const &first, 
        typename token_sequence_type::const_iterator const &last, 
        typename token_type::position_type const &tok,
        bool if_block_status);
};

///////////////////////////////////////////////////////////////////////////////
}   //  namespace grammars
}   //  namespace wave 
}   //  namespace boost

#endif // !defined(CPP_EXPRESSION_GRAMMAR_GEN_HPP_42399258_6CDC_4101_863D_5C7D95B5A6CA_INCLUDED)
