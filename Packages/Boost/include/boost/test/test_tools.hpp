//  (C) Copyright Gennadiy Rozental 2001-2003.
//  (C) Copyright Ullrich Koethe 2001.
//  Use, modification, and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : contains definition for all test tools in test toolbox
// ***************************************************************************

#ifndef BOOST_TEST_TOOLS_HPP
#define BOOST_TEST_TOOLS_HPP

// Boost.Test
#include <boost/test/detail/unit_test_config.hpp>
#include <boost/test/unit_test_log.hpp>
#include <boost/test/detail/class_properties.hpp>
#include <boost/test/detail/wrap_stringstream.hpp>

// BOOST
#include <boost/cstdlib.hpp> // for boost::exit_success;
#include <boost/config.hpp>  // compilers workarounds
#include <boost/shared_ptr.hpp>

#include <stdexcept>        // for std::exception
#include <cstddef>          // for std::size_t
#include <memory>           // for std::auto_ptr
#include <string>           // for std::string

// ************************************************************************** //
// **************                    TOOL BOX                  ************** //
// ************************************************************************** //

#define BOOST_CHECKPOINT(message_)                                          \
    boost::test_toolbox::detail::checkpoint_impl(                           \
        boost::wrap_stringstream().ref() << message_, __FILE__, __LINE__)   \
/**/

#define BOOST_WARN(predicate)                                               \
    boost::test_toolbox::detail::warn_and_continue_impl((predicate),        \
        boost::wrap_stringstream().ref() << #predicate, __FILE__, __LINE__) \
/**/

#define BOOST_CHECK(predicate)                                              \
    boost::test_toolbox::detail::test_and_continue_impl((predicate),        \
        boost::wrap_stringstream().ref() << #predicate, __FILE__, __LINE__) \
/**/

#define BOOST_CHECK_EQUAL(left_, right_)                                                \
    boost::test_toolbox::detail::equal_and_continue_impl((left_), (right_),             \
        boost::wrap_stringstream().ref() << #left_ " == " #right_, __FILE__, __LINE__)
/**/

#define BOOST_CHECK_CLOSE(left_, right_, tolerance) \
    boost::test_toolbox::detail::compare_and_continue_impl((left_), (right_), (tolerance),\
                                                           #left_,  #right_, __FILE__, __LINE__)
/**/

#define BOOST_BITWISE_EQUAL(left_, right_) \
    boost::test_toolbox::detail::bitwise_equal_and_continue_impl((left_), (right_), \
        boost::wrap_stringstream().ref() << #left_ " =.= " #right_, __FILE__, __LINE__)
/**/

#define BOOST_REQUIRE(predicate) \
    boost::test_toolbox::detail::test_and_throw_impl((predicate), \
        boost::wrap_stringstream().ref() << #predicate, __FILE__, __LINE__)
/**/

#define BOOST_MESSAGE(message_) \
    boost::test_toolbox::detail::message_impl( \
        boost::wrap_stringstream().ref() << message_, __FILE__, __LINE__)
/**/

#define BOOST_WARN_MESSAGE(predicate, message_) \
    boost::test_toolbox::detail::warn_and_continue_impl((predicate), \
        boost::wrap_stringstream().ref() << message_, __FILE__, __LINE__,false)
/**/

#define BOOST_CHECK_MESSAGE(predicate, message_) \
    boost::test_toolbox::detail::test_and_continue_impl((predicate), \
        boost::wrap_stringstream().ref() << message_, __FILE__, __LINE__,false)
/**/

#define BOOST_REQUIRE_MESSAGE(predicate, message_) \
    boost::test_toolbox::detail::test_and_throw_impl((predicate), \
        boost::wrap_stringstream().ref() << message_, __FILE__, __LINE__,false)
/**/

#define BOOST_CHECK_PREDICATE( predicate, arg_list_size, arg_list ) \
    boost::test_toolbox::detail::test_and_continue_impl(predicate, BOOST_PLACE_PREDICATE_ARGS ## arg_list_size arg_list, \
        boost::wrap_stringstream().ref() << #predicate << "("\
        << BOOST_PRINT_PREDICATE_ARGS ## arg_list_size arg_list << ")", __FILE__, __LINE__)
/**/

#define BOOST_REQUIRE_PREDICATE( predicate, arg_list_size, arg_list ) \
    boost::test_toolbox::detail::test_and_throw_impl(predicate, BOOST_PLACE_PREDICATE_ARGS ## arg_list_size arg_list, \
        boost::wrap_stringstream().ref() << #predicate << "("\
        << BOOST_PRINT_PREDICATE_ARGS ## arg_list_size arg_list << ")", __FILE__, __LINE__)
/**/

#define BOOST_ERROR(message_) BOOST_CHECK_MESSAGE( false, message_ )

#define BOOST_FAIL(message_) BOOST_REQUIRE_MESSAGE( false, message_ )

#define BOOST_CHECK_THROW( statement, exception )                                               \
    try { statement; BOOST_ERROR( "exception "#exception" is expected" ); }                     \
    catch( exception const& ) {                                                                 \
        BOOST_CHECK_MESSAGE( true, "exception "#exception" is caught" );                        \
    }                                                                                           \
/**/

#define BOOST_CHECK_EXCEPTION( statement, exception, predicate )                                \
    try { statement; BOOST_ERROR( "exception "#exception" is expected" ); }                     \
    catch( exception const& ex ) {                                                              \
        BOOST_CHECK_MESSAGE( predicate( ex ), "incorrect exception "#exception" is caught" );   \
    }                                                                                           \
/**/

#define BOOST_IGNORE_CHECK( e ) true

#define BOOST_CHECK_NO_THROW( statement ) \
    try { statement; BOOST_CHECK_MESSAGE( true, "no exceptions was thrown by "#statement ); } \
    catch( ... ) { \
        BOOST_ERROR( "exception was thrown by "#statement ); \
    }
/**/

#define BOOST_CHECK_EQUAL_COLLECTIONS( left_begin_, left_end_, right_begin_ ) \
    boost::test_toolbox::detail::equal_and_continue_impl( (left_begin_), (left_end_), (right_begin_),\
        boost::wrap_stringstream().ref() << \
            "{" #left_begin_ ", " #left_end_ "}" " == {" #right_begin_ ", ...}", __FILE__, __LINE__)
/**/

#define BOOST_IS_DEFINED(symb) boost::test_toolbox::detail::is_defined_impl( #symb, BOOST_STRINGIZE(= symb) )

// ***************************** //
// helper macros

#define BOOST_PLACE_PREDICATE_ARGS1( first_ ) first_
#define BOOST_PLACE_PREDICATE_ARGS2( first_, second_ ) first_, second_

#define BOOST_PRINT_PREDICATE_ARGS1( first_ ) #first_
#define BOOST_PRINT_PREDICATE_ARGS2( first_, second_ ) #first_ << ", " << #second_

// ***************************** //
// depricated interface

#define BOOST_TEST(predicate)            BOOST_CHECK(predicate)
#define BOOST_CRITICAL_TEST(predicate)   BOOST_REQUIRE(predicate)
#define BOOST_CRITICAL_ERROR(message_)   BOOST_FAIL(message_)

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4511) // copy constructor could not be generated
# pragma warning(disable: 4512) // assignment operator could not be generated
#endif

namespace boost {

namespace test_toolbox {

// ************************************************************************** //
// **************            extended_predicate_value          ************** //
// ************************************************************************** //

class extended_predicate_value {
public:
    // Constructor
    extended_predicate_value( bool predicate_value_ )
    : p_predicate_value( predicate_value_ ), p_message( new wrap_stringstream ) {}

    extended_predicate_value( extended_predicate_value const& rhs )
    : p_predicate_value( rhs.p_predicate_value.get() ), 
      p_message( rhs.p_message )                    {}

    bool        operator!() const                   { return !p_predicate_value.get(); }
    void        operator=( bool predicate_value_ )  { p_predicate_value.value = predicate_value_; }

    BOOST_READONLY_PROPERTY( bool, 1, (extended_predicate_value) )
                p_predicate_value;
    boost::shared_ptr<wrap_stringstream>
                p_message;
};

namespace detail {

using unit_test_framework::c_string_literal;

// ************************************************************************** //
// **************                test_tool_failed              ************** //
// ************************************************************************** //

// exception used to implement critical checks

struct test_tool_failed : public std::exception {
};

// ************************************************************************** //
// **************               log print helper               ************** //
// ************************************************************************** //

template<typename T>
struct print_log_value {
    void    operator()( std::ostream& ostr, T const& t )
    {
        ostr << t; // by default print the value
    }
};

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

#define BOOST_TEST_DONT_PRINT_LOG_VALUE( the_type )                 \
namespace boost { namespace test_toolbox { namespace detail {       \
template<>                                                          \
struct print_log_value<the_type > {                                 \
    void operator()( std::ostream& ostr, the_type const& t ) {}     \
};                                                                  \
}}}                                                                 \
/**/

//____________________________________________________________________________//

template<typename T>
struct print_helper {
    explicit    print_helper( T const& t ) : m_t( t ) {}

    T const&    m_t;
};

//____________________________________________________________________________//

template<typename T>
inline std::ostream& 
operator<<( std::ostream& ostr, print_helper<T> const& ph )
{
    print_log_value<T>()( ostr, ph.m_t );

    return ostr;
}

//____________________________________________________________________________//

// ************************************************************************** //
// **************            TOOL BOX Implementation           ************** //
// ************************************************************************** //

void
checkpoint_impl( wrap_stringstream& message_, c_string_literal file_name_, std::size_t line_num_ );

//____________________________________________________________________________//

void
message_impl( wrap_stringstream& message_, c_string_literal file_name_, std::size_t line_num_ );

//____________________________________________________________________________//

// ************************************* //

void
warn_and_continue_impl( bool predicate_, wrap_stringstream& message_,
                        c_string_literal file_name_, std::size_t line_num_,
                        bool add_fail_pass_ = true );

//____________________________________________________________________________//

void
warn_and_continue_impl( extended_predicate_value const& v_, wrap_stringstream& message_,
                        c_string_literal file_name_, std::size_t line_num_,
                        bool add_fail_pass_ = true );

//____________________________________________________________________________//

// ************************************* //

bool  // return true if error detected
test_and_continue_impl( bool predicate_, wrap_stringstream& message_,
                        c_string_literal file_name_, std::size_t line_num_,
                        bool add_fail_pass_ = true,
                        unit_test_framework::log_level log_level_ = unit_test_framework::log_all_errors );
void
test_and_throw_impl   ( bool predicate_, wrap_stringstream& message_,
                        c_string_literal file_name_, std::size_t line_num_,
                        bool add_fail_pass_ = true,
                        unit_test_framework::log_level log_level_ = unit_test_framework::log_fatal_errors );

//____________________________________________________________________________//

bool
test_and_continue_impl( extended_predicate_value const& v_, wrap_stringstream& message_,
                        c_string_literal file_name_, std::size_t line_num_,
                        bool add_fail_pass_ = true,
                        unit_test_framework::log_level log_level_ = unit_test_framework::log_all_errors );

//____________________________________________________________________________//

// Borland bug workaround
#if BOOST_WORKAROUND(__BORLANDC__, <= 0x570)
inline bool
test_and_continue_impl( void* ptr, wrap_stringstream& message_,
                        c_string_literal file_name_, std::size_t line_num_,
                        bool add_fail_pass_ = true,
                        unit_test_framework::log_level log_level_ = unit_test_framework::log_all_errors )
{
    return test_and_continue_impl( !!ptr, message_, file_name_, line_num_, add_fail_pass_, log_level_ );
}
#endif

//____________________________________________________________________________//

void
test_and_throw_impl   ( extended_predicate_value const& v_, wrap_stringstream& message_,
                        c_string_literal file_name_, std::size_t line_num_,
                        bool add_fail_pass_ = true,
                        unit_test_framework::log_level log_level_ = unit_test_framework::log_fatal_errors );

//____________________________________________________________________________//

template<typename ArgType, typename Predicate>
inline bool
test_and_continue_impl( Predicate const& pred_, ArgType const& arg_,
                        wrap_stringstream& message_,
                        c_string_literal file_name_, std::size_t line_num_,
                        unit_test_framework::log_level log_level_ = unit_test_framework::log_all_errors )
{
    extended_predicate_value predicate( pred_( arg_ ) );

    if( !predicate ) {
        return test_and_continue_impl( predicate,
                                       wrap_stringstream().ref() << "test " << message_ << " failed for " 
                                                                 << print_helper<ArgType>( arg_ ),
                                       file_name_, line_num_, false, log_level_ );
    }

    return test_and_continue_impl( predicate, message_, file_name_, line_num_, true, log_level_ );
}

//____________________________________________________________________________//

template<typename ArgType, typename Predicate>
inline void
test_and_throw_impl   ( Predicate const& pred_, ArgType const& arg_,
                        wrap_stringstream& message_,
                        c_string_literal file_name_, std::size_t line_num_,
                        unit_test_framework::log_level log_level_ = unit_test_framework::log_fatal_errors )
{
    if( test_and_continue_impl( arg_, pred_, message_, file_name_, line_num_, log_level_ ) ) {
        throw test_tool_failed(); // error already reported by test_and_continue_impl
    }
}

//____________________________________________________________________________//

template<typename First, typename Second, typename Predicate>
inline bool
test_and_continue_impl( Predicate const& pred_, First const& first_, Second const& second_,
                        wrap_stringstream& message_,
                        c_string_literal file_name_, std::size_t line_num_,
                        unit_test_framework::log_level log_level_ = unit_test_framework::log_all_errors )
{
    extended_predicate_value predicate( pred_( first_, second_ ) );

    if( !predicate ) {
        return test_and_continue_impl( predicate,
            wrap_stringstream().ref() << "test " << message_ 
                                      << " failed for (" << print_helper<First>( first_ ) << ", " 
                                                         << print_helper<Second>( second_ ) << ")",
            file_name_, line_num_, false, log_level_ );
    }

    return test_and_continue_impl( predicate, message_, file_name_, line_num_, true, log_level_ );
}

//____________________________________________________________________________//

template<typename First, typename Second, typename Predicate>
inline void
test_and_throw_impl( First const& first_, Second const& second_, Predicate const& pred_,
                     wrap_stringstream& message_, c_string_literal file_name_, std::size_t line_num_,
                     unit_test_framework::log_level log_level_ = unit_test_framework::log_fatal_errors )
{
    if( test_and_continue_impl( first_, second_, pred_, message_, file_name_, line_num_, log_level_ ) ) {
        throw test_tool_failed(); // error already reported by test_and_continue_impl
    }
}

//____________________________________________________________________________//

// ************************************* //

bool
equal_and_continue_impl( c_string_literal left_, c_string_literal right_, wrap_stringstream& message_,
                         c_string_literal file_name_, std::size_t line_num_,
                         unit_test_framework::log_level log_level_ = unit_test_framework::log_all_errors );

//____________________________________________________________________________//

template <class Left, class Right>
inline bool
equal_and_continue_impl( Left const& left_, Right const& right_,
                         wrap_stringstream& message_, c_string_literal file_name_, std::size_t line_num_,
                         unit_test_framework::log_level log_level_ = unit_test_framework::log_all_errors,
                         std::size_t pos = (std::size_t)-1 )
{
    extended_predicate_value predicate( left_ == right_ );

    if( !predicate ) {
        wrap_stringstream error_message;
        error_message.ref() << "test " << message_ << " failed";

        if( pos != (std::size_t)-1 )
            error_message.ref() <<  " in a position " << pos;

        error_message.ref() << " [" 
                            << print_helper<Left>( left_ )   << " != " 
                            << print_helper<Right>( right_ ) << "]";

        return test_and_continue_impl( predicate, error_message, file_name_, line_num_, false, log_level_ );
    }

    return test_and_continue_impl( predicate, wrap_stringstream().ref() << message_, file_name_, line_num_, true, log_level_ );
    //----------------------------------------------^ this is added to prevent message_ corruption when reused by collection comparison
}

//____________________________________________________________________________//

template <class Left, class Right>
inline void
equal_and_continue_impl( Left left_begin_, Left left_end_, Right right_begin_,
                         wrap_stringstream& message_,
                         c_string_literal file_name_, std::size_t line_num_,
                         unit_test_framework::log_level log_level_ = unit_test_framework::log_all_errors )
{
    std::size_t pos = 0;
    for( ; left_begin_ != left_end_; ++left_begin_, ++right_begin_, ++pos )
        equal_and_continue_impl( *left_begin_, *right_begin_, message_, file_name_, line_num_, log_level_, pos );
}

//____________________________________________________________________________//

// ************************************* //

template<typename FPT, typename PersentType>
inline bool
compare_and_continue_impl( FPT left_, FPT right_, PersentType tolerance_,
                           c_string_literal left_text_, c_string_literal right_text_,
                           c_string_literal file_name_, std::size_t line_num_,
                           unit_test_framework::log_level log_level_ = unit_test_framework::log_all_errors )
{
    extended_predicate_value predicate( check_is_closed( left_, right_, tolerance_ ) );

    if( !predicate ) {
        return test_and_continue_impl( predicate,
            wrap_stringstream().ref() << "difference between " << left_text_ << "{" << print_helper<FPT>( left_ ) << "}" 
                                      << " and " << right_text_ << "{" << print_helper<FPT>( right_ ) << "}" 
                                      << " exceeds " << print_helper<PersentType>( tolerance_ ) << "%",
            file_name_, line_num_, false, log_level_ );
    }

    return test_and_continue_impl( predicate, 
        wrap_stringstream().ref() << "difference between " << left_text_ << "{" << print_helper<FPT>( left_ ) << "}" 
                                  << " and " << right_text_ << "{" << print_helper<FPT>( right_ ) << "}" 
                                  << " does not exceeds " << print_helper<PersentType>( tolerance_ ) << "%",
        file_name_, line_num_, true, log_level_ );
}

//____________________________________________________________________________//

template <class Left, class Right>
inline void
bitwise_equal_and_continue_impl( Left const& left_, Right const& right_,
                                 wrap_stringstream& message_, char const* file_name_, std::size_t line_num_,
                                 unit_test_framework::log_level log_level_ = unit_test_framework::log_all_errors )
{
    std::size_t left_bit_size  = sizeof(Left)*CHAR_BIT;
    std::size_t right_bit_size = sizeof(Right)*CHAR_BIT;

    static Left const L1( 1 );
    static Right const R1( 1 );

    if( left_bit_size != right_bit_size )
        warn_and_continue_impl( false, wrap_stringstream().ref() << message_ << ": operands bit sizes does not coinside", 
                                file_name_, line_num_, false );

    std::size_t total_bits = left_bit_size < right_bit_size ? left_bit_size : right_bit_size;

    for( std::size_t counter = 0; counter < total_bits; ++counter ) {
        bool predicate = ( left_ & ( L1 << counter ) ) == ( right_ & ( R1 << counter ) );

        test_and_continue_impl( predicate, wrap_stringstream().ref() << message_.str() << " in the position " << counter,
                                file_name_, line_num_, true, log_level_ );
    }
}

//____________________________________________________________________________//

// ************************************* //

bool
is_defined_impl( c_string_literal symbol_name_, c_string_literal symbol_value_ );

//____________________________________________________________________________//

} // namespace detail

// ************************************************************************** //
// **************               output_test_stream             ************** //
// ************************************************************************** //

// class to be used to simplify testing of ostream print functions

class output_test_stream : public 
#ifdef BOOST_NO_STRINGSTREAM
    std::ostrstream
#else
    std::ostringstream
#endif // BOOST_NO_STRINGSTREAM
{
    typedef extended_predicate_value                result_type;
    typedef unit_test_framework::c_string_literal   c_string_literal;
public:
    // Constructor
    explicit        output_test_stream( std::string const&  pattern_file_name = std::string(),
                                        bool                match_or_save     = true );
    explicit        output_test_stream( c_string_literal    pattern_file_name,
                                        bool                match_or_save     = true );

    // Destructor
    ~output_test_stream();

    // checking function
    result_type     is_empty( bool flush_stream_ = true );
    result_type     check_length( std::size_t length_, bool flush_stream_ = true );
    result_type     is_equal( c_string_literal arg_, bool flush_stream_ = true );
    result_type     is_equal( std::string const& arg_, bool flush_stream_ = true );
    result_type     is_equal( c_string_literal arg_, std::size_t n_, bool flush_stream_ = true );
    result_type     match_pattern( bool flush_stream_ = true );

    // helper function
    void            flush();
    std::size_t     length();

private:
    void            sync();

    struct Impl;
    boost::shared_ptr<Impl> m_pimpl;
};

} // namespace test_toolbox

} // namespace boost

#ifdef BOOST_MSVC
# pragma warning(default: 4511) // copy constructor could not be generated
# pragma warning(default: 4512) // assignment operator could not be generated
# pragma warning(pop)
#endif

// ***************************************************************************
//  Revision History :
//
//  $Log$
//  Revision 1.1  2004/05/23 10:51:36  spion
//  Initial revision
//
//  Revision 1.35.2.1  2004/01/06 13:33:28  johnmaddock
//  merged changes from main branch
//
//  Revision 1.36  2004/01/05 11:56:25  johnmaddock
//  Borland specific workaround needs to be inline to prevent linker errors, and is unneeded for version 6 of the compiler.
//
//  Revision 1.35  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//

// ***************************************************************************

#endif // BOOST_TEST_TOOLS_HPP
