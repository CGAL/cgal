//  Boost integer/static_log2.hpp header file  -------------------------------//

//  (C) Copyright Daryle Walker 2001.  Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears 
//  in all copies.  This software is provided "as is" without express or
//  implied warranty, and with no claim as to its suitability for any purpose. 

//  See http://www.boost.org for updates, documentation, and revision history. 

#ifndef BOOST_INTEGER_STATIC_LOG2_HPP
#define BOOST_INTEGER_STATIC_LOG2_HPP

#include <boost/integer_fwd.hpp>  // self include

#include <boost/config.hpp>  // for BOOST_STATIC_CONSTANT, etc.
#include <boost/limits.hpp>  // for std::numeric_limits

#ifdef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
#include <boost/pending/ct_if.hpp>  // for boost::ct_if<>
#endif


namespace boost
{


//  Implementation details  --------------------------------------------------//

namespace detail
{

// Forward declarations
template < unsigned long Val, int Place = 0, int Index
 = std::numeric_limits<unsigned long>::digits >
    struct static_log2_helper_t;

#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

template < unsigned long Val, int Place >
    struct static_log2_helper_t< Val, Place, 1 >;

#else

template < int Place >
    struct static_log2_helper_final_step;

template < unsigned long Val, int Place = 0, int Index
 = std::numeric_limits<unsigned long>::digits >
    struct static_log2_helper_nopts_t;

#endif

// Recursively build the logarithm by examining the upper bits
template < unsigned long Val, int Place, int Index >
struct static_log2_helper_t
{
private:
    BOOST_STATIC_CONSTANT( int, half_place = Index / 2 );
    BOOST_STATIC_CONSTANT( unsigned long, lower_mask = (1ul << half_place)
     - 1ul );
    BOOST_STATIC_CONSTANT( unsigned long, upper_mask = ~lower_mask );
    BOOST_STATIC_CONSTANT( bool, do_shift = (Val & upper_mask) != 0ul );

    BOOST_STATIC_CONSTANT( unsigned long, new_val = do_shift ? (Val
     >> half_place) : Val );
    BOOST_STATIC_CONSTANT( int, new_place = do_shift ? (Place + half_place)
     : Place );
    BOOST_STATIC_CONSTANT( int, new_index = Index - half_place );

#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
    typedef static_log2_helper_t<new_val, new_place, new_index>  next_step_type;
#else
    typedef static_log2_helper_nopts_t<new_val, new_place, new_index>  next_step_type;
#endif

public:
    BOOST_STATIC_CONSTANT( int, value = next_step_type::value );

};  // boost::detail::static_log2_helper_t

// Non-recursive case
#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

template < unsigned long Val, int Place >
struct static_log2_helper_t< Val, Place, 1 >
{
public:
    BOOST_STATIC_CONSTANT( int, value = Place );

};  // boost::detail::static_log2_helper_t

#else

template < int Place >
struct static_log2_helper_final_step
{
public:
    BOOST_STATIC_CONSTANT( int, value = Place );

};  // boost::detail::static_log2_helper_final_step

template < unsigned long Val, int Place, int Index >
struct static_log2_helper_nopts_t
{
private:
    typedef static_log2_helper_t<Val, Place, Index>  recursive_step_type;
    typedef static_log2_helper_final_step<Place>     final_step_type;

    typedef typename ct_if<( Index != 1 ), recursive_step_type,
     final_step_type>::type  next_step_type;

public:
    BOOST_STATIC_CONSTANT( int, value = next_step_type::value );

};  // boost::detail::static_log2_helper_nopts_t

#endif

}  // namespace detail


//  Compile-time log-base-2 evaluator class declaration  ---------------------//

template < unsigned long Value >
struct static_log2
{
    BOOST_STATIC_CONSTANT( int, value
     = detail::static_log2_helper_t<Value>::value );
};

template < >
struct static_log2< 0ul >
{
    // The logarithm of zero is undefined.
};


}  // namespace boost


#endif  // BOOST_INTEGER_STATIC_LOG2_HPP
