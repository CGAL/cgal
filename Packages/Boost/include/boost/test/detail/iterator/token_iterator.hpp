//  (C) Copyright Gennadiy Rozental 2004.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : token iterator for string and range tokenization
// ***************************************************************************

#ifndef BOOST_TOKEN_ITERATOR_HPP_071894GER
#define BOOST_TOKEN_ITERATOR_HPP_071894GER

// Boost
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/iterator_traits.hpp>

#include <boost/test/detail/iterator/input_iterator_facade.hpp>
#include <boost/test/detail/basic_cstring/basic_cstring.hpp>

// STL
#include <iosfwd>
#include <cctype>

#ifdef BOOST_NO_STDC_NAMESPACE
namespace std{ using ::ispunct; using ::isspace; }
#endif

namespace boost {

namespace unit_test {

// ************************************************************************** //
// **************               ti_delimeter_type              ************** //
// ************************************************************************** //

enum ti_delimeter_type { use_delim, use_ispunct, use_isspace };

namespace ut_detail {

// ************************************************************************** //
// **************                 delim_policy                 ************** //
// ************************************************************************** //

template<typename CharT,typename CharCompare>
class delim_policy {
    typedef basic_cstring<CharT const> cstring;
public:
    // Constructor
    explicit    delim_policy( ti_delimeter_type t = use_delim ) : m_type( t ) {}
    explicit    delim_policy( cstring d, ti_delimeter_type t )
    {
        use_delimeters( d, t );
    }

    void        use_delimeters( cstring d, ti_delimeter_type t )
    {
        m_delimeters = d;
        m_type       = d.is_empty() ? t : use_delim;
    }

    bool        operator()( CharT c )
    {
        switch( m_type ) {
        case use_delim: {
            typename cstring::iterator it = m_delimeters.begin();
            for( ; it != m_delimeters.end(); ++it )
                if( CharCompare()( *it, c ) )
                    break;
            return it != m_delimeters.end();
        }
        case use_ispunct:   return (std::ispunct)( c ) != 0;
        case use_isspace:   return (std::isspace)( c ) != 0;
        }

        return false;
    }

    // Data members
    cstring             m_delimeters;
    ti_delimeter_type   m_type;
};

// ************************************************************************** //
// **************                 token_assigner               ************** //
// ************************************************************************** //

template<typename TraversalTag>
struct token_assigner {
#if BOOST_WORKAROUND( BOOST_DINKUMWARE_STDLIB, < 306 )
    template<typename Iterator, typename C, typename T>
    static void assign( Iterator b, Iterator e, std::basic_string<C,T>& t )
    { for( ; b != e; ++b ) t += *b; }
    
    template<typename Iterator, typename C>
    static void assign( Iterator b, Iterator e, basic_cstring<C>& t )  { t.assign( b, e ); }
#else
    template<typename Iterator, typename Token>
    static void assign( Iterator b, Iterator e, Token& t )  { t.assign( b, e ); }
#endif
    template<typename Iterator, typename Token>
    static void append_move( Iterator& b, Token& )          { ++b; }

    template<typename Token>
    static void clear( Token& )                             {}
};

template<>
struct token_assigner<single_pass_traversal_tag> {
    template<typename Iterator, typename Token>
    static void assign( Iterator b, Iterator e, Token& t )  {}

    template<typename Iterator, typename Token>
    static void append_move( Iterator& b, Token& t )        { t += *b; ++b; }

#if BOOST_WORKAROUND( BOOST_DINKUMWARE_STDLIB, < 306 ) \
    || BOOST_WORKAROUND( __GNUC__, < 3 ) && !defined(__SGI_STL_PORT) \
        && !defined(_STLPORT_VERSION)
    template<typename Token>
    static void clear( Token& t )                           { t.erase(); }
#else
    template<typename Token>
    static void clear( Token& t )                           { t.clear(); }
#endif
};

// ************************************************************************** //
// **************             default_char_compare             ************** //
// ************************************************************************** //

template<typename CharT>
class default_char_compare {
public:
    bool operator()( CharT c1, CharT c2 )
    {
#ifdef BOOST_CLASSIC_IOSTREAMS
        return std::string_char_traits<CharT>::eq( c1, c2 );
#else
        return std::char_traits<CharT>::eq( c1, c2 );
#endif
    }
};

// ************************************************************************** //
// **************             dropped_delimeters_m             ************** //
// ************************************************************************** //

template<typename CharT>
struct dropped_delimeters_m
{
    explicit dropped_delimeters_m( basic_cstring<CharT const> delims )
    : m_delims( delims ) {}

    template<typename CharCompare>
    void apply( delim_policy<CharT,CharCompare>& is_drop, delim_policy<CharT,CharCompare>&, bool ) const
    {
        is_drop.use_delimeters( m_delims, use_isspace );
    }

    basic_cstring<CharT const> m_delims;
};

//____________________________________________________________________________//

template<>
struct dropped_delimeters_m<ti_delimeter_type>
{
    explicit dropped_delimeters_m( ti_delimeter_type type )
    : m_type( type ) {}

    template<typename CharT, typename CharCompare>
    void apply( delim_policy<CharT,CharCompare>& is_drop, delim_policy<CharT,CharCompare>&, bool ) const
    {
        is_drop.use_delimeters( "", m_type );
    }

    ti_delimeter_type m_type;
};

//____________________________________________________________________________//

// ************************************************************************** //
// **************               kept_delimeters_m              ************** //
// ************************************************************************** //

template<typename CharT>
struct kept_delimeters_m
{
    explicit kept_delimeters_m( basic_cstring<CharT const> delims, ti_delimeter_type type = use_ispunct )
    : m_delims( delims ), m_type( type ) {}

    template<typename CharCompare>
    void apply( delim_policy<CharT,CharCompare>&, delim_policy<CharT,CharCompare>& is_kept, bool ) const
    {
        is_kept.use_delimeters( m_delims, m_type );
    }

    basic_cstring<CharT const>  m_delims;
    ti_delimeter_type           m_type;
};

//____________________________________________________________________________//

template<>
struct kept_delimeters_m<ti_delimeter_type>
{
    explicit kept_delimeters_m( ti_delimeter_type type )
    : m_type( type ) {}

    template<typename CharT, typename CharCompare>
    void apply( delim_policy<CharT,CharCompare>&, delim_policy<CharT,CharCompare>& is_kept, bool ) const
    {
        is_kept.use_delimeters( "", m_type );
    }

    ti_delimeter_type m_type;
};

//____________________________________________________________________________//

// ************************************************************************** //
// **************              keep_empty_tokens_m             ************** //
// ************************************************************************** //

struct keep_empty_tokens_m
{
    explicit keep_empty_tokens_m( bool v )
    : m_keep_empty_tokens( v ) {}

    template<typename CharT,typename CharCompare>
    void apply( delim_policy<CharT,CharCompare>&, delim_policy<CharT,CharCompare>&, bool& keep_empty_tokens ) const
    {
        keep_empty_tokens = m_keep_empty_tokens;
    }

    bool m_keep_empty_tokens;
};

} // namespace ut_detail

// ************************************************************************** //
// **************              modifiers generators            ************** //
// ************************************************************************** //

static struct dropped_delimeters_generator {
    template<typename CharT>
    ut_detail::dropped_delimeters_m<CharT>
    operator=( basic_cstring<CharT const> d ) { return ut_detail::dropped_delimeters_m<CharT>( d ); }
    template<typename CharT>
    ut_detail::dropped_delimeters_m<CharT>
    operator=( CharT const* d ) { return ut_detail::dropped_delimeters_m<CharT>( basic_cstring<CharT const>( d ) ); }

    ut_detail::dropped_delimeters_m<ti_delimeter_type>
    operator=( ti_delimeter_type t ) { return ut_detail::dropped_delimeters_m<ti_delimeter_type>( t ); }
} dropped_delimeters;

//____________________________________________________________________________//

static struct kept_delimeters_generator {
    template<typename CharT>
    ut_detail::kept_delimeters_m<CharT>
    operator=( basic_cstring<CharT const> d ) { return ut_detail::kept_delimeters_m<CharT>( d ); }
    template<typename CharT>
    ut_detail::kept_delimeters_m<CharT>
    operator=( CharT const* d ) { return ut_detail::kept_delimeters_m<CharT>( basic_cstring<CharT const>( d ) ); }

    ut_detail::kept_delimeters_m<ti_delimeter_type>
    operator=( ti_delimeter_type t ) { return ut_detail::kept_delimeters_m<ti_delimeter_type>( t ); }
} kept_delimeters;

//____________________________________________________________________________//

static struct keep_empty_tokens_generator : ut_detail::keep_empty_tokens_m {
    keep_empty_tokens_generator() : ut_detail::keep_empty_tokens_m( true ) {}

    ut_detail::keep_empty_tokens_m
    operator=( bool v ) { return ut_detail::keep_empty_tokens_m( v ); }
} keep_empty_tokens;

//____________________________________________________________________________//

// ************************************************************************** //
// **************             token_iterator_base              ************** //
// ************************************************************************** //

template<typename Derived,
         typename CharT,
         typename CharCompare   = ut_detail::default_char_compare<CharT>,
         typename ValueType     = basic_cstring<CharT const>,
         typename Reference     = basic_cstring<CharT const>,
         typename Traversal     = forward_traversal_tag>
class token_iterator_base
: public input_iterator_facade<Derived,ValueType,Reference,Traversal> {
    typedef basic_cstring<CharT const>                                      cstring;
    typedef ut_detail::delim_policy<CharT,CharCompare>                      delim_policy;
    typedef input_iterator_facade<Derived,ValueType,Reference,Traversal>    base;

public:
    // Constructor
    explicit    token_iterator_base()
    : m_is_dropped( use_isspace ),
      m_is_kept( use_ispunct ),
      m_keep_empty_tokens( false ),
      m_token_produced( false )
    {
    }

protected:
    template<typename Modifier>
    token_iterator_base&
    use( Modifier const& m )
    {
        m.apply( m_is_dropped, m_is_kept, m_keep_empty_tokens );
        return *this;
    }

    template<typename Iter> 
    bool                    get( Iter& begin, Iter end )
    {
        typedef ut_detail::token_assigner<BOOST_DEDUCED_TYPENAME iterator_traversal<Iter>::type> Assigner;
        Iter checkpoint;

        Assigner::clear( this->m_value );

        if( !m_keep_empty_tokens ) {
            while( begin != end && m_is_dropped( *begin ) )
                ++begin;

            if( begin == end )
                return false;

            checkpoint = begin;

            if( m_is_kept( *begin ) )
                Assigner::append_move( begin, this->m_value );
            else
                while( begin != end && !m_is_dropped( *begin ) && !m_is_kept( *begin ) )
                    Assigner::append_move( begin, this->m_value );
        } 
        else { // m_keep_empty_tokens ia true
            checkpoint = begin;

            if( begin == end ) {
                if( m_token_produced ) 
                    return false;

                m_token_produced = true;
            }
            else if( m_is_kept( *begin ) ) {
                if( m_token_produced ) 
                    Assigner::append_move( begin, this->m_value );

                m_token_produced = !m_token_produced;
            } 
            else if( !m_token_produced && m_is_dropped( *begin ) )
                m_token_produced = true;
            else {
                if( m_is_dropped( *begin ) )
                    checkpoint = ++begin;

                while( begin != end && !m_is_dropped( *begin ) && !m_is_kept( *begin ) )
                    Assigner::append_move( begin, this->m_value );

                m_token_produced = true;
            }
        }

        Assigner::assign( checkpoint, begin, this->m_value );

        return true;
    }

private:
    // Data members
    delim_policy            m_is_dropped;
    delim_policy            m_is_kept;
    bool                    m_keep_empty_tokens;
    bool                    m_token_produced;
};

// ************************************************************************** //
// **************          basic_string_token_iterator         ************** //
// ************************************************************************** //

template<typename CharT,
         typename CharCompare = ut_detail::default_char_compare<CharT> >
class basic_string_token_iterator
: public token_iterator_base<basic_string_token_iterator<CharT,CharCompare>,CharT,CharCompare> {
    typedef basic_cstring<CharT const> cstring;
    typedef token_iterator_base<basic_string_token_iterator<CharT,CharCompare>,CharT,CharCompare> base;
public:
    explicit    basic_string_token_iterator() {}
    explicit    basic_string_token_iterator( cstring src )
    : m_src( src )
    {
        this->init();
    }

    template<typename Src, typename M1>
    basic_string_token_iterator( Src src, M1 const& m1 )
    : m_src( src )
    {
        this->use( m1 );

        this->init();
    }

    template<typename Src, typename M1, typename M2>
    basic_string_token_iterator( Src src, M1 const& m1, M2 const& m2 )
    : m_src( src )
    {
        this->use( m1 );
        this->use( m2 );

        this->init();
    }

    template<typename Src, typename M1, typename M2, typename M3>
    basic_string_token_iterator( Src src, M1 const& m1, M2 const& m2, M3 const& m3 )
    : m_src( src )
    {
        this->use( m1 );
        this->use( m2 );
        this->use( m3 );

        this->init();
    }

private:
    friend class input_iterator_core_access;

    // input iterator implementation
    bool        get()
    {
        typename cstring::iterator begin = m_src.begin();
        bool res = base::get( begin, m_src.end() );

        m_src.assign( begin, m_src.end() );

        return res;
    }

    // Data members
    cstring     m_src;
};

typedef basic_string_token_iterator<char>       string_token_iterator;
typedef basic_string_token_iterator<wchar_t>    wstring_token_iterator;

// ************************************************************************** //
// **************              range_token_iterator            ************** //
// ************************************************************************** //

template<typename Iter,
         typename CharCompare = ut_detail::default_char_compare<BOOST_DEDUCED_TYPENAME iterator_value<Iter>::type>,
         typename ValueType   = std::basic_string<BOOST_DEDUCED_TYPENAME iterator_value<Iter>::type>,
         typename Reference   = ValueType const&>
class range_token_iterator
: public token_iterator_base<range_token_iterator<Iter,CharCompare,ValueType,Reference>,
                             typename iterator_value<Iter>::type,CharCompare,ValueType,Reference> {
    typedef basic_cstring<typename ValueType::value_type> cstring;
    typedef token_iterator_base<range_token_iterator<Iter,CharCompare,ValueType,Reference>,
                                typename iterator_value<Iter>::type,CharCompare,ValueType,Reference> base;
public:
    explicit    range_token_iterator() {}
    explicit    range_token_iterator( Iter begin, Iter end = Iter() )
    : m_begin( begin ), m_end( end )
    {
        this->init();
    }

    template<typename M1>
    range_token_iterator( Iter begin, Iter end, M1 const& m1 )
    : m_begin( begin ), m_end( end )
    {
        this->use( m1 );

        this->init();
    }

    template<typename M1,typename M2>
    range_token_iterator( Iter begin, Iter end, M1 const& m1, M2 const& m2 )
    : m_begin( begin ), m_end( end )
    {
        this->use( m1 );
        this->use( m2 );

        this->init();
    }

    template<typename M1,typename M2,typename M3>
    range_token_iterator( Iter begin, Iter end, M1 const& m1, M2 const& m2, M3 const& m3 )
    : m_begin( begin ), m_end( end )
    {
        this->use( m1 );
        this->use( m2 );
        this->use( m3 );

        this->init();
    }

private:
    friend class input_iterator_core_access;

    // input iterator implementation
    bool        get()
    {
        return base::get( m_begin, m_end );
    }

    // Data members
    Iter m_begin;
    Iter m_end;
};

// ************************************************************************** //
// **************            make_range_token_iterator         ************** //
// ************************************************************************** //

template<typename Iter>
inline range_token_iterator<Iter>
make_range_token_iterator( Iter begin, Iter end = Iter() )
{
    return range_token_iterator<Iter>( begin, end );
}

//____________________________________________________________________________//

template<typename Iter,typename M1>
inline range_token_iterator<Iter>
make_range_token_iterator( Iter begin, Iter end, M1 const& m1 )
{
    return range_token_iterator<Iter>( begin, end, m1 );
}

//____________________________________________________________________________//

template<typename Iter, typename M1, typename M2>
inline range_token_iterator<Iter>
make_range_token_iterator( Iter begin, Iter end, M1 const& m1, M2 const& m2 )
{
    return range_token_iterator<Iter>( begin, end, m1, m2 );
}

//____________________________________________________________________________//

template<typename Iter,typename M1, typename M2, typename M3>
inline range_token_iterator<Iter>
make_range_token_iterator( Iter begin, Iter end, M1 const& m1, M2 const& m2, M3 const& m3 )
{
    return range_token_iterator<Iter>( begin, end, m1, m2, m3 );
}

//____________________________________________________________________________//

} // namespace unit_test

} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/11/20 10:52:25  spion
//  Initial revision
//
//  Revision 1.7.2.1  2004/10/30 11:34:27  agurtovoy
//  MSVC/Borland fixes
//
//  Revision 1.7  2004/10/01 10:50:40  rogeeff
//  gcc 2.95 workarounds
//
//  Revision 1.6  2004/09/27 08:38:08  rogeeff
//  msvc workaround
//
//  Revision 1.5  2004/09/19 09:22:13  rogeeff
//  ios fix for classic iostreams
//
//  Revision 1.4  2004/07/20 12:46:39  vladimir_prus
//  Add "this->", since gcc 3.4 has two-phase lookup.
//
//  Revision 1.3  2004/07/19 12:29:57  rogeeff
//  guard rename
//  mingw port
//
//  Revision 1.2  2004/06/07 07:33:50  rogeeff
//  detail namespace renamed
//
//  Revision 1.1  2004/06/05 11:03:12  rogeeff
//  input_iterator_adaptor simplified
//  token_iterator added
//
// ***************************************************************************

#endif // BOOST_TOKEN_ITERATOR_HPP_071894GER

