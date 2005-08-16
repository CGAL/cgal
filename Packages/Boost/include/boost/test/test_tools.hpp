//  (C) Copyright Gennadiy Rozental 2005.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : contains definition for all test tools in test toolbox
// ***************************************************************************

#ifndef BOOST_TEST_TEST_TOOLS_HPP_012705GER
#define BOOST_TEST_TEST_TOOLS_HPP_012705GER

// Boost.Test
#include <boost/test/predicate_result.hpp>

#include <boost/test/detail/config.hpp>
#include <boost/test/detail/global_typedef.hpp>
#include <boost/test/detail/workaround.hpp>

#include <boost/test/utils/wrap_stringstream.hpp>
#include <boost/test/utils/basic_cstring/io.hpp>

// Boost
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/seq/enum.hpp> 
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/limits.hpp>

#include <boost/type_traits/is_array.hpp>
#include <boost/type_traits/is_function.hpp>

#include <boost/mpl/or.hpp>

// STL
#include <cstddef>          // for std::size_t
#include <iosfwd>

#include <boost/test/detail/suppress_warnings.hpp>

//____________________________________________________________________________//

// ************************************************************************** //
// **************                    TOOL BOX                  ************** //
// ************************************************************************** //

// In macros below following argument abbreviations are used:
// P - predicate
// M - message
// S - statement
// E - exception
// L - left argument
// R - right argument
// TL - tool level
// CT - check type
// ARGS - arguments list

#define BOOST_TEST_TOOL_IMPL( func, P, check_descr, TL, CT ) \
    boost::test_tools::tt_detail::func(                      \
        P,                                                   \
        boost::wrap_stringstream().ref() << check_descr,     \
        BOOST_TEST_L(__FILE__),                              \
        (std::size_t)__LINE__,                               \
        boost::test_tools::tt_detail::TL,                    \
        boost::test_tools::tt_detail::CT                     \
/**/

//____________________________________________________________________________//

#define BOOST_CHECK_IMPL( P, check_descr, TL, CT )  BOOST_TEST_TOOL_IMPL( check_impl, P, check_descr, TL, CT ), 0 )

//____________________________________________________________________________//

#define BOOST_TEST_PASS_ARG_INFO( r, data, arg ) , arg, BOOST_STRINGIZE( arg )

#define BOOST_CHECK_WITH_ARGS_IMPL( P, check_descr, TL, CT, ARGS )  \
    BOOST_TEST_TOOL_IMPL( check_frwd, P, check_descr, TL, CT )      \
    BOOST_PP_SEQ_FOR_EACH( BOOST_TEST_PASS_ARG_INFO, '_', ARGS ) )  \
/**/

//____________________________________________________________________________//

#define BOOST_WARN( P )                     BOOST_CHECK_IMPL( (P), BOOST_TEST_STRINGIZE( P ), WARN, CHECK_PRED )
#define BOOST_CHECK( P )                    BOOST_CHECK_IMPL( (P), BOOST_TEST_STRINGIZE( P ), CHECK, CHECK_PRED )
#define BOOST_REQUIRE( P )                  BOOST_CHECK_IMPL( (P), BOOST_TEST_STRINGIZE( P ), REQUIRE, CHECK_PRED )

//____________________________________________________________________________//

#define BOOST_WARN_MESSAGE( P, M )          BOOST_CHECK_IMPL( (P), M, WARN, CHECK_MSG )
#define BOOST_CHECK_MESSAGE( P, M )         BOOST_CHECK_IMPL( (P), M, CHECK, CHECK_MSG )
#define BOOST_REQUIRE_MESSAGE( P, M )       BOOST_CHECK_IMPL( (P), M, REQUIRE, CHECK_MSG )

//____________________________________________________________________________//

#define BOOST_ERROR( M )                    BOOST_CHECK_MESSAGE( false, M )
#define BOOST_FAIL( M )                     BOOST_REQUIRE_MESSAGE( false, M )

//____________________________________________________________________________//

#define BOOST_MESSAGE( M )                  BOOST_CHECK_IMPL( false, M, WARN, MSG_ONLY )

//____________________________________________________________________________//

#define BOOST_CHECKPOINT( M )               BOOST_CHECK_IMPL( false, M, WARN, SET_CHECKPOINT )

//____________________________________________________________________________//

#define BOOST_CHECK_THROW_IMPL( S, E, P, prefix, TL )                                                   \
    try {                                                                                               \
        S;                                                                                              \
        BOOST_CHECK_IMPL( false, "exception " BOOST_STRINGIZE( E ) " is expected", TL, CHECK_MSG ); }   \
    catch( E const& ex ) {                                                                              \
        boost::unit_test::ut_detail::ignore_unused_variable_warning( ex );                              \
        BOOST_CHECK_IMPL( P, prefix BOOST_STRINGIZE( E ) " is caught", TL, CHECK_MSG );                 \
    }                                                                                                   \
/**/

//____________________________________________________________________________//

#define BOOST_WARN_THROW( S, E )            BOOST_CHECK_THROW_IMPL( S, E, true, "exception ", WARN )
#define BOOST_CHECK_THROW( S, E )           BOOST_CHECK_THROW_IMPL( S, E, true, "exception ", CHECK )
#define BOOST_REQUIRE_THROW( S, E )         BOOST_CHECK_THROW_IMPL( S, E, true, "exception ", REQUIRE )

//____________________________________________________________________________//

#define BOOST_WARN_EXCEPTION( S, E, P )     BOOST_CHECK_THROW_IMPL( S, E, P( ex ), "incorrect exception ", WARN )
#define BOOST_CHECK_EXCEPTION( S, E, P )    BOOST_CHECK_THROW_IMPL( S, E, P( ex ), "incorrect exception ", CHECK )
#define BOOST_REQUIRE_EXCEPTION( S, E, P )  BOOST_CHECK_THROW_IMPL( S, E, P( ex ), "incorrect exception ", REQUIRE )

//____________________________________________________________________________//

#define BOOST_IGNORE_CHECK( e )             true

//____________________________________________________________________________//

#define BOOST_CHECK_NO_THROW_IMPL( S, TL )                                                          \
    try {                                                                                           \
        S;                                                                                          \
        BOOST_CHECK_IMPL( true, "no exceptions thrown by " BOOST_STRINGIZE( S ), TL, CHECK_MSG ); } \
    catch( ... ) {                                                                                  \
        BOOST_CHECK_IMPL( false, "exception thrown by " BOOST_STRINGIZE( S ), TL, CHECK_MSG );      \
    }                                                                                               \
/**/

#define BOOST_WARN_NO_THROW( S )            BOOST_CHECK_NO_THROW_IMPL( S, WARN )
#define BOOST_CHECK_NO_THROW( S )           BOOST_CHECK_NO_THROW_IMPL( S, CHECK )
#define BOOST_REQUIRE_NO_THROW( S )         BOOST_CHECK_NO_THROW_IMPL( S, REQUIRE )

//____________________________________________________________________________//

#define BOOST_WARN_EQUAL( L, R ) \
    BOOST_CHECK_WITH_ARGS_IMPL( boost::test_tools::tt_detail::equal_impl_frwd(), "", WARN, CHECK_EQUAL, (L)(R) )
#define BOOST_CHECK_EQUAL( L, R ) \
    BOOST_CHECK_WITH_ARGS_IMPL( boost::test_tools::tt_detail::equal_impl_frwd(), "", CHECK, CHECK_EQUAL, (L)(R) )
#define BOOST_REQUIRE_EQUAL( L, R ) \
    BOOST_CHECK_WITH_ARGS_IMPL( boost::test_tools::tt_detail::equal_impl_frwd(), "", REQUIRE, CHECK_EQUAL, (L)(R) )

//____________________________________________________________________________//

#define BOOST_WARN_CLOSE( L, R, T ) \
    BOOST_CHECK_WITH_ARGS_IMPL( boost::test_tools::check_is_close, "", WARN, CHECK_CLOSE, (L)(R)(T) )
#define BOOST_CHECK_CLOSE( L, R, T ) \
    BOOST_CHECK_WITH_ARGS_IMPL( boost::test_tools::check_is_close, "", CHECK, CHECK_CLOSE, (L)(R)(T) )
#define BOOST_REQUIRE_CLOSE( L, R, T ) \
    BOOST_CHECK_WITH_ARGS_IMPL( boost::test_tools::check_is_close, "", REQUIRE, CHECK_CLOSE, (L)(R)(T) )

//____________________________________________________________________________//

#define BOOST_WARN_SMALL( FPV, T ) \
    BOOST_CHECK_WITH_ARGS_IMPL( boost::test_tools::check_is_small, "", WARN, CHECK_SMALL, (FPV)(T) )
#define BOOST_CHECK_SMALL( FPV, T ) \
    BOOST_CHECK_WITH_ARGS_IMPL( boost::test_tools::check_is_small, "", CHECK, CHECK_SMALL, (FPV)(T) )
#define BOOST_REQUIRE_SMALL( FPV, T ) \
    BOOST_CHECK_WITH_ARGS_IMPL( boost::test_tools::check_is_small, "", REQUIRE, CHECK_SMALL, (FPV)(T) )

//____________________________________________________________________________//

#define BOOST_WARN_PREDICATE( P, ARGS ) \
    BOOST_CHECK_WITH_ARGS_IMPL( P, BOOST_TEST_STRINGIZE( P ), WARN, CHECK_PRED_WITH_ARGS, ARGS )
#define BOOST_CHECK_PREDICATE( P, ARGS ) \
    BOOST_CHECK_WITH_ARGS_IMPL( P, BOOST_TEST_STRINGIZE( P ), CHECK, CHECK_PRED_WITH_ARGS, ARGS )
#define BOOST_REQUIRE_PREDICATE( P, ARGS ) \
    BOOST_CHECK_WITH_ARGS_IMPL( P, BOOST_TEST_STRINGIZE( P ), REQUIRE, CHECK_PRED_WITH_ARGS, ARGS )

//____________________________________________________________________________//

#define BOOST_EQUAL_COLLECTIONS_IMPL( L_begin, L_end, R_begin, R_end, TL )      \
    BOOST_TEST_TOOL_IMPL( check_impl, boost::test_tools::tt_detail::equal_coll_impl( \
        (L_begin), (L_end), (R_begin), (R_end) ), "", TL, CHECK_EQUAL_COLL ),   \
    4,                                                                          \
    BOOST_STRINGIZE( L_begin ), BOOST_STRINGIZE( L_end ),                       \
    BOOST_STRINGIZE( R_begin ), BOOST_STRINGIZE( R_end ) )                      \
/**/

#define BOOST_WARN_EQUAL_COLLECTIONS( L_begin, L_end, R_begin, R_end )          \
    BOOST_EQUAL_COLLECTIONS_IMPL( L_begin, L_end, R_begin, R_end, WARN )
#define BOOST_CHECK_EQUAL_COLLECTIONS( L_begin, L_end, R_begin, R_end )         \
    BOOST_EQUAL_COLLECTIONS_IMPL( L_begin, L_end, R_begin, R_end, CHECK )
#define BOOST_REQUIRE_EQUAL_COLLECTIONS( L_begin, L_end, R_begin, R_end )       \
    BOOST_EQUAL_COLLECTIONS_IMPL( L_begin, L_end, R_begin, R_end, REQUIRE )

//____________________________________________________________________________//

#define BOOST_BITWISE_EQUAL_IMPL( L, R, TL )                                    \
    BOOST_TEST_TOOL_IMPL( check_impl,                                           \
      boost::test_tools::tt_detail::bitwise_equal_impl( (L), (R) ),             \
      "", TL, CHECK_BITWISE_EQUAL ),                                            \
    2, BOOST_STRINGIZE( L ), BOOST_STRINGIZE( R ) )                             \
/**/

#define BOOST_WARN_BITWISE_EQUAL( L, R )    BOOST_BITWISE_EQUAL_IMPL( L, R, WARN )
#define BOOST_CHECK_BITWISE_EQUAL( L, R )   BOOST_BITWISE_EQUAL_IMPL( L, R, CHECK )
#define BOOST_REQUIRE_BITWISE_EQUAL( L, R ) BOOST_BITWISE_EQUAL_IMPL( L, R, REQUIRE )

//____________________________________________________________________________//

#define BOOST_IS_DEFINED( symb )            boost::test_tools::tt_detail::is_defined_impl( #symb, BOOST_STRINGIZE(= symb) )

//____________________________________________________________________________//

// ***************************** //
// deprecated interface

#define BOOST_BITWISE_EQUAL( L, R )         BOOST_CHECK_BITWISE_EQUAL( L, R )

namespace boost {

namespace test_tools {

typedef unit_test::const_string      const_string;

namespace tt_detail {

// ************************************************************************** //
// **************              tools classification            ************** //
// ************************************************************************** //

enum check_type {
    CHECK_PRED, 
    CHECK_MSG,
    CHECK_EQUAL,
    CHECK_CLOSE,
    CHECK_SMALL,
    CHECK_BITWISE_EQUAL,
    MSG_ONLY,
    SET_CHECKPOINT, 
    CHECK_PRED_WITH_ARGS,
    CHECK_EQUAL_COLL
};

enum tool_level {
    WARN, CHECK, REQUIRE, PASS
};

// ************************************************************************** //
// **************               log print helper               ************** //
// ************************************************************************** //

template<typename T>
struct print_log_value {
    void    operator()( std::ostream& ostr, T const& t )
    {
        typedef typename mpl::or_<is_array<T>,is_function<T> >::type couldnt_use_nl;

        set_precision( ostr, couldnt_use_nl() );

        ostr << t; // by default print the value
    }

    void set_precision( std::ostream& ostr, mpl::false_ )
    {
        if( std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::radix == 2 )
            ostr.precision( 2 + std::numeric_limits<T>::digits * 301/1000 ); 
    }

    void set_precision( std::ostream&, mpl::true_ ) {}
};

//____________________________________________________________________________//

#define BOOST_TEST_DONT_PRINT_LOG_VALUE( the_type )                 \
namespace boost { namespace test_tools { namespace tt_detail {      \
template<>                                                          \
struct print_log_value<the_type > {                                 \
    void operator()( std::ostream& ostr, the_type const& t ) {}     \
};                                                                  \
}}}                                                                 \
/**/

//____________________________________________________________________________//

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
template<typename T, std::size_t N >
struct print_log_value< T[N] > {
    void    operator()( std::ostream& ostr, T const* t )
    {   
        ostr << t;
    }
};
#endif

//____________________________________________________________________________//

template<>
struct print_log_value<char> {
    void    operator()( std::ostream& ostr, char t );
};

//____________________________________________________________________________//

template<>
struct print_log_value<unsigned char> {
    void    operator()( std::ostream& ostr, unsigned char t );
};

//____________________________________________________________________________//

template<>
struct print_log_value<char const*> {
    void    operator()( std::ostream& ostr, char const* t );
};

//____________________________________________________________________________//

template<>
struct print_log_value<wchar_t const*> {
    void    operator()( std::ostream& ostr, wchar_t const* t );
};

//____________________________________________________________________________//

template<typename T>
struct print_helper_t {
    explicit    print_helper_t( T const& t ) : m_t( t ) {}

    T const&    m_t;
};

//____________________________________________________________________________//

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
// Borland suffers premature pointer decay passing arrays by reference
template<typename T, std::size_t N >
struct print_helper_t< T[N] > {
    explicit    print_helper_t( T const * t ) : m_t( t ) {}

    T const *   m_t;
};
#endif

//____________________________________________________________________________//

template<typename T>
inline print_helper_t<T> print_helper( T const& t )
{
    return print_helper_t<T>( t );
}

//____________________________________________________________________________//

template<typename T>
inline std::ostream& 
operator<<( std::ostream& ostr, print_helper_t<T> const& ph )
{
    print_log_value<T>()( ostr, ph.m_t );

    return ostr;
}

//____________________________________________________________________________//

// ************************************************************************** //
// **************            TOOL BOX Implementation           ************** //
// ************************************************************************** //

void check_impl( predicate_result const& pr, wrap_stringstream& check_descr,
                 const_string file_name, std::size_t line_num,
                 tool_level tool_level, check_type ct,
                 std::size_t num_args, ... );

//____________________________________________________________________________//

#define TEMPL_PARAMS( z, m, dummy ) , typename BOOST_JOIN( Arg, m )
#define FUNC_PARAMS( z, m, dummy )                                                  \
    , BOOST_JOIN( Arg, m ) const& BOOST_JOIN( arg, m )                              \
    , char const* BOOST_JOIN( BOOST_JOIN( arg, m ), _descr )                        \
/**/

#define PRED_PARAMS( z, m, dummy ) BOOST_PP_COMMA_IF( m ) BOOST_JOIN( arg, m ) 

#define ARG_INFO( z, m, dummy )                                                     \
    , BOOST_JOIN( BOOST_JOIN( arg, m ), _descr )                                    \
    , (boost::wrap_stringstream().ref()                                             \
        << boost::test_tools::tt_detail::                                           \
            print_helper( BOOST_JOIN( arg, m ) )).str().c_str()                     \
/**/

#define IMPL_FRWD( z, n, dummy )                                                    \
template<typename Pred                                                              \
         BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), TEMPL_PARAMS, _ )>            \
inline void                                                                         \
check_frwd( Pred P, wrap_stringstream& check_descr,                                 \
            const_string file_name, std::size_t line_num,                           \
            tool_level tool_level, check_type ct                                    \
            BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), FUNC_PARAMS, _ )           \
)                                                                                   \
{                                                                                   \
    check_impl( P( BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), PRED_PARAMS, _ ) ), \
                check_descr, file_name, line_num, tool_level, ct,                   \
                BOOST_PP_ADD( n, 1 )                                                \
                BOOST_PP_REPEAT_ ## z( BOOST_PP_ADD( n, 1 ), ARG_INFO, _ )          \
    );                                                                              \
}                                                                                   \
/**/

#ifndef BOOST_TEST_MAX_PREDICATE_ARITY
#define BOOST_TEST_MAX_PREDICATE_ARITY 5
#endif

BOOST_PP_REPEAT( BOOST_TEST_MAX_PREDICATE_ARITY, IMPL_FRWD, _ )

#undef TEMPL_PARAMS
#undef FUNC_PARAMS
#undef PRED_INFO
#undef ARG_INFO
#undef IMPL_FRWD

//____________________________________________________________________________//

template <class Left, class Right>
predicate_result    equal_impl( Left const& left, Right const& right )
{
    return left == right;
}

//____________________________________________________________________________//

predicate_result        equal_impl( char const* left, char const* right );
inline predicate_result equal_impl( char* left, char const* right ) { return equal_impl( (char const*)left, (char const*)right ); }
inline predicate_result equal_impl( char const* left, char* right ) { return equal_impl( (char const*)left, (char const*)right ); }
inline predicate_result equal_impl( char* left, char* right )       { return equal_impl( (char const*)left, (char const*)right ); }

#if !defined( BOOST_NO_CWCHAR )
predicate_result        equal_impl( wchar_t const* left, wchar_t const* right );
inline predicate_result equal_impl( wchar_t* left, wchar_t const* right ) { return equal_impl( (wchar_t const*)left, (wchar_t const*)right ); }
inline predicate_result equal_impl( wchar_t const* left, wchar_t* right ) { return equal_impl( (wchar_t const*)left, (wchar_t const*)right ); }
inline predicate_result equal_impl( wchar_t* left, wchar_t* right )       { return equal_impl( (wchar_t const*)left, (wchar_t const*)right ); }
#endif

//____________________________________________________________________________//

struct equal_impl_frwd {
    template <typename Left, typename Right>
    inline predicate_result
    call_impl( Left const& left, Right const& right, mpl::false_ ) const
    {
        return equal_impl( left, right );
    }

    template <typename Left, typename Right>
    inline predicate_result
    call_impl( Left const& left, Right const& right, mpl::true_ ) const
    {
        return (*this)( right, &left[0] );
    }

    template <typename Left, typename Right>
    inline predicate_result
    operator()( Left const& left, Right const& right ) const
    {
        typedef typename is_array<Left>::type left_is_array;
        return call_impl( left, right, left_is_array() );
    }
};

//____________________________________________________________________________//

template <typename Left, typename Right>
inline predicate_result
equal_coll_impl( Left left_begin, Left left_end, Right right_begin, Right right_end )
{
    predicate_result    res( true );
    std::size_t         pos = 0;

    for( ; left_begin != left_end && right_begin != right_end; ++left_begin, ++right_begin, ++pos ) {
        if( *left_begin != *right_begin ) {
            res = false;
            res.message() << "\nMismatch in a position " << pos << ": "  << *left_begin << " != " << *right_begin;
        }
    }

    if( left_begin != left_end ) {
        std::size_t r_size = pos;
        while( left_begin != left_end ) {
            ++pos;
            ++left_begin;
        }

        res = false;
        res.message() << "\nCollections size mismatch: " << pos << " != " << r_size;
    }

    if( right_begin != right_end ) {
        std::size_t l_size = pos;
        while( right_begin != right_end ) {
            ++pos;
            ++right_begin;
        }

        res = false;
        res.message() << "\nCollections size mismatch: " << l_size << " != " << pos;
    }

    return res;
}

//____________________________________________________________________________//

template <class Left, class Right>
inline predicate_result
bitwise_equal_impl( Left const& left, Right const& right )
{
    predicate_result    res( true );

    std::size_t left_bit_size  = sizeof(Left)*CHAR_BIT;
    std::size_t right_bit_size = sizeof(Right)*CHAR_BIT;

    static Left const  L1( 1 );
    static Right const R1( 1 );

    std::size_t total_bits = left_bit_size < right_bit_size ? left_bit_size : right_bit_size;

    for( std::size_t counter = 0; counter < total_bits; ++counter ) {
        if( ( left & ( L1 << counter ) ) != ( right & ( R1 << counter ) ) ) {
            res = false;
            res.message() << "\nMismatch in a position " << counter;
        }
    }

    if( left_bit_size != right_bit_size ) {
        res = false;
        res.message() << "\nOperands bit sizes mismatch: " << left_bit_size << " != " << right_bit_size;
    }

    return res;
}

//____________________________________________________________________________//

bool is_defined_impl( const_string symbol_name, const_string symbol_value );

//____________________________________________________________________________//

} // namespace tt_detail

} // namespace test_tools

namespace test_toolbox = test_tools;

} // namespace boost

//____________________________________________________________________________//

#include <boost/test/detail/enable_warnings.hpp>

// ***************************************************************************
//  Revision History :
//
//  $Log$
//  Revision 1.1.1.3  2005/08/16 11:24:10  spion
//  Import of Boost v. 1.33.0
//
//  Revision 1.54  2005/06/07 04:38:20  rogeeff
//  borland fix
//
//  Revision 1.53  2005/05/11 04:51:14  rogeeff
//  borlard portability fix
//
//  Revision 1.52  2005/03/22 07:08:47  rogeeff
//  string comparisons streamlined
//  precision settings made portable
//
//  Revision 1.51  2005/02/21 10:23:54  rogeeff
//  major issue with TT redesign causing TT to reevaluate it's arguments fixed
//  FP precision extended
//
//  Revision 1.50  2005/02/20 08:27:06  rogeeff
//  This a major update for Boost.Test framework. See release docs for complete list of fixes/updates
//
//  Revision 1.49  2005/02/01 06:40:06  rogeeff
//  copyright update
//  old log entries removed
//  minor stylistic changes
//  deprecated tools removed
//
//  Revision 1.48  2005/01/30 03:32:57  rogeeff
//  Test Tools completely reworked:
//    interfaces streamlined to provide 3 version for each tool
//    implementation reworked to use single vararg formatter function
//    CHECK_COLLECTION now expect 4 arguments
//    BITWISE_EQUAL renamed to CHECK_BITWISE_EQUAL but still provided as deprecated
//    CHECK_COLLECTION interface changed to use PP_SEQ and as a result support arbitrary number of predicate arguments
//    most of templates eliminated
//    deprecated tools removed
//    print_helper object generator added
//
// ***************************************************************************

#endif // BOOST_TEST_TEST_TOOLS_HPP_012705GER
