/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    A generic C++ lexer token definition
    
    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(CPP_TOKEN_HPP_53A13BD2_FBAA_444B_9B8B_FCB225C2BBA8_INCLUDED)
#define CPP_TOKEN_HPP_53A13BD2_FBAA_444B_9B8B_FCB225C2BBA8_INCLUDED

#include <boost/wave/wave_config.hpp>
#include <boost/wave/util/file_position.hpp>
#include <boost/wave/token_ids.hpp>  
#include <boost/wave/language_support.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace cpplexer {

///////////////////////////////////////////////////////////////////////////////
//  forward declaration of the token type
template <typename PositionT = boost::wave::util::file_position_type>
class lex_token;

///////////////////////////////////////////////////////////////////////////////
//  
//  lex_token
//
///////////////////////////////////////////////////////////////////////////////

template <typename PositionT>
class lex_token 
{
public:
    typedef BOOST_WAVE_STRINGTYPE   string_type;
    typedef PositionT               position_type;
    
    lex_token()
    :   id(T_EOI)
    {}
    
    lex_token(token_id id_, string_type const &value_, PositionT const &pos_)
    :   id(id_), value(value_), pos(pos_)
    {}

// accessors
    operator token_id() const { return id; }
    string_type const &get_value() const { return value; }
    position_type const &get_position() const { return pos; }
    void set_token_id (token_id id_) { id = id_; }
    void set_value (string_type const &newval) { value = newval; }
    void set_position (position_type const &pos_) { pos = pos_; }

    friend bool operator== (lex_token const& lhs, lex_token const& rhs)
    {
        //  two tokens are considered equal even if they contain different 
        //  positions
        return (lhs.id == rhs.id && lhs.value == rhs.value) ? true : false;
    }
    
// debug support    
#if BOOST_WAVE_DUMP_PARSE_TREE != 0
// access functions for the tree_to_xml functionality
    static int get_token_id(lex_token const &t) 
        { return ID_FROM_TOKEN(token_id(t)); }
    static string_type get_token_value(lex_token const &t) 
        { return t.get_value(); }
#endif 
    
#if defined(BOOST_SPIRIT_DEBUG)
// debug support
    void print (std::ostream &stream) const
    {
        stream << get_token_name(id) << "(";
        for (std::size_t i = 0; i < value.size(); ++i) {
            switch (value[i]) {
            case '\r':  stream << "\\r"; break;
            case '\n':  stream << "\\n"; break;
            default:
                stream << value[i]; 
                break;
            }
        }
        stream << ")";
    }
#endif // defined(BOOST_SPIRIT_DEBUG)

private:
    token_id id;                // the token id
    string_type value;          // the text, which was parsed into this token
    PositionT pos;              // the original file position
};

#if defined(BOOST_SPIRIT_DEBUG)
template <typename PositionT>
inline std::ostream &
operator<< (std::ostream &stream, lex_token<PositionT> const &object)
{
    object.print(stream);
    return stream;
}
#endif // defined(BOOST_SPIRIT_DEBUG)

///////////////////////////////////////////////////////////////////////////////
}   // namespace cpplexer
}   // namespace wave
}   // namespace boost

#endif // !defined(CPP_TOKEN_HPP_53A13BD2_FBAA_444B_9B8B_FCB225C2BBA8_INCLUDED)
