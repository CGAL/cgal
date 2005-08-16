/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library
    Whitespace eater
    
    http://www.boost.org/

    Copyright (c) 2003 Paul Mensonides
    Copyright (c) 2001-2005 Hartmut Kaiser. 
    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(EAT_WHITESPACE_HPP_4CE9AD17_F82D_4AB2_A117_555DF0DCC801_INCLUDED)
#define EAT_WHITESPACE_HPP_4CE9AD17_F82D_4AB2_A117_555DF0DCC801_INCLUDED

#include <boost/wave/wave_config.hpp>   
#include <boost/wave/token_ids.hpp>   

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace util {

///////////////////////////////////////////////////////////////////////////////
template <typename TokenT>
bool ccomment_has_newline(TokenT const& token)
{
    using namespace boost::wave;

    if (T_CCOMMENT == token_id(token)) {
        if (TokenT::string_type::npos != 
            token.get_value().find_first_of("\n"))
        {
            return true;
        }
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////////
template <typename TokenT>
class eat_whitespace {

public:
    eat_whitespace(bool preserve_comments);
    
    bool may_skip (TokenT &token, bool &skipped_newline);

protected:
    bool skip_cppcomment(boost::wave::token_id id)
    {
        return !preserve_comments && T_CPPCOMMENT == id;
    }
    
private:
    typedef bool state_t(TokenT &token, bool &skipped_newline);
    state_t eat_whitespace::* state;
    state_t general, newline, newline_2nd, whitespace;
    bool preserve_comments;
};

template <typename TokenT>
inline 
eat_whitespace<TokenT>::eat_whitespace(bool preserve_comments)
:   state(&eat_whitespace::newline), preserve_comments(preserve_comments)
{
}

template <typename TokenT>
inline bool 
eat_whitespace<TokenT>::may_skip(TokenT &token, bool &skipped_newline) 
{
    return (this->*state)(token, skipped_newline);
}

template <typename TokenT>
inline bool 
eat_whitespace<TokenT>::general(TokenT &token, bool &skipped_newline) 
{
    using namespace boost::wave;

    token_id id = token_id(token);
    if (T_NEWLINE == id || T_CPPCOMMENT == id) {
        state = &eat_whitespace::newline;
    }
    else if (T_SPACE == id || T_SPACE2 == id || T_CCOMMENT == id) {
        state = &eat_whitespace::whitespace;

        if (ccomment_has_newline(token)) 
            skipped_newline = true;

        if ((!preserve_comments || T_CCOMMENT != id) && 
            token.get_value().size() > 1)
        {
            token.set_value(" ");   // replace with a single space
        }
    }
    else {
        state = &eat_whitespace::general;
    }
    return false;
}

template <typename TokenT>
inline bool 
eat_whitespace<TokenT>::newline(TokenT &token, bool &skipped_newline) 
{
    using namespace boost::wave;
    
    token_id id = token_id(token);
    if (T_NEWLINE == id || T_CPPCOMMENT == id) {
        skipped_newline = true;
        state = &eat_whitespace::newline_2nd;
        return T_NEWLINE == id || skip_cppcomment(id);
    }
    else if (T_SPACE != id && T_SPACE2 != id && T_CCOMMENT != id) {
        return general(token, skipped_newline);
    }

    if (T_CCOMMENT == id) {
        if (ccomment_has_newline(token))
            skipped_newline = true;

        if (preserve_comments) {
            state = &eat_whitespace::general;
            return false;
        }
        // fall through...
    }
    return true;
}

template <typename TokenT>
inline bool 
eat_whitespace<TokenT>::newline_2nd(TokenT &token, bool &skipped_newline) 
{
    using namespace boost::wave;
    
    token_id id = token_id(token);
    if (T_SPACE == id || T_SPACE2 == id)
        return true;
    if (T_CCOMMENT == id) {
        if (ccomment_has_newline(token))
            skipped_newline = true;

        if (preserve_comments) {
            state = &eat_whitespace::general;
            return false;
        }
        return  true;
    }
    if (T_NEWLINE != id && T_CPPCOMMENT != id) 
        return general(token, skipped_newline);

    skipped_newline = true;
    return T_NEWLINE == id || skip_cppcomment(id);
}

template <typename TokenT>
inline bool 
eat_whitespace<TokenT>::whitespace(TokenT &token, bool &skipped_newline) 
{
    using namespace boost::wave;
    
    token_id id = token_id(token);
    if (T_SPACE != id && T_SPACE2 != id && 
        T_CCOMMENT != id && T_CPPCOMMENT != id) 
    {
        return general(token, skipped_newline);
    }
    
    if (T_CCOMMENT == id) {
        if (ccomment_has_newline(token))
            skipped_newline = true;
        return !preserve_comments;
    }

    return T_SPACE == id || T_SPACE2 == id || skip_cppcomment(id);
}

///////////////////////////////////////////////////////////////////////////////
}   // namespace util
}   // namespace wave
}   // namespace boost

#endif // !defined(EAT_WHITESPACE_HPP_4CE9AD17_F82D_4AB2_A117_555DF0DCC801_INCLUDED)

