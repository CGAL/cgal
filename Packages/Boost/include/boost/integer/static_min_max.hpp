//  Boost integer/static_min_max.hpp header file  ----------------------------//

//  (C) Copyright Daryle Walker 2001.  Permission to copy, use, modify, sell
//  and distribute this software is granted provided this copyright notice
//  appears in all copies.  This software is provided "as is" without
//  express or implied warranty, and with no claim as to its suitability
//  for any purpose. 

//  See http://www.boost.org for updates, documentation, and revision history. 

#ifndef BOOST_INTEGER_STATIC_MIN_MAX_HPP
#define BOOST_INTEGER_STATIC_MIN_MAX_HPP

#include <boost/integer_fwd.hpp>  // self include

#include <boost/config.hpp>  // for BOOST_STATIC_CONSTANT


namespace boost
{


//  Compile-time extrema class declarations  ---------------------------------//
//  Get the minimum or maximum of two values, signed or unsigned.

template < long Value1, long Value2 >
struct static_signed_min
{
    BOOST_STATIC_CONSTANT( long, value = (Value1 > Value2) ? Value2 : Value1 );
};

template < long Value1, long Value2 >
struct static_signed_max
{
    BOOST_STATIC_CONSTANT( long, value = (Value1 < Value2) ? Value2 : Value1 );
};

template < unsigned long Value1, unsigned long Value2 >
struct static_unsigned_min
{
    BOOST_STATIC_CONSTANT( unsigned long, value
     = (Value1 > Value2) ? Value2 : Value1 );
};

template < unsigned long Value1, unsigned long Value2 >
struct static_unsigned_max
{
    BOOST_STATIC_CONSTANT( unsigned long, value
     = (Value1 < Value2) ? Value2 : Value1 );
};


}  // namespace boost


#endif  // BOOST_INTEGER_STATIC_MIN_MAX_HPP
