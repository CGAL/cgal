//  boost utility/base_from_member.hpp header file  --------------------------//

//  Copyright 2001, 2003 Daryle Walker.  Use, modification, and distribution are
//  subject to the Boost Software License, Version 1.0.  (See accompanying file
//  LICENSE_1_0.txt or a copy at <http://www.boost.org/LICENSE_1_0.txt>.)

//  See <http://www.boost.org/libs/utility/> for the library's home page.

#ifndef BOOST_UTILITY_BASE_FROM_MEMBER_HPP
#define BOOST_UTILITY_BASE_FROM_MEMBER_HPP

#include <boost/utility_fwd.hpp>  // required for parameter defaults


namespace boost
{

//  Base-from-member class template  -----------------------------------------//

// Helper to initialize a base object so a derived class can use this
// object in the initialization of another base class.  Used by
// Dietmar Kuehl from ideas by Ron Klatcho to solve the problem of a
// base class needing to be initialized by a member.

// Contributed by Daryle Walker

template < typename MemberType, int UniqueID >
class base_from_member
{
protected:
    MemberType  member;

    base_from_member()
        : member()
        {}

    template< typename T1 >
    explicit base_from_member( T1 x1 )
        : member( x1 )
        {}

    template< typename T1, typename T2 >
    base_from_member( T1 x1, T2 x2 )
        : member( x1, x2 )
        {}

    template< typename T1, typename T2, typename T3 >
    base_from_member( T1 x1, T2 x2, T3 x3 )
        : member( x1, x2, x3 ) 
        {}

    template< typename T1, typename T2, typename T3, typename T4 >
    base_from_member( T1 x1, T2 x2, T3 x3, T4 x4 )
        : member( x1, x2, x3, x4 ) 
        {}

    template< typename T1, typename T2, typename T3, typename T4, typename T5 >
    base_from_member( T1 x1, T2 x2, T3 x3, T4 x4, T5 x5 )
        : member( x1, x2, x3, x4, x5 ) 
        {}

    template< typename T1, typename T2, typename T3, typename T4, typename T5,
     typename T6 >
    base_from_member( T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6 )
        : member( x1, x2, x3, x4, x5, x6 ) 
        {}

    template< typename T1, typename T2, typename T3, typename T4, typename T5,
     typename T6, typename T7 >
    base_from_member( T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7 )
        : member( x1, x2, x3, x4, x5, x6, x7 ) 
        {}

    template< typename T1, typename T2, typename T3, typename T4, typename T5,
     typename T6, typename T7, typename T8 >
    base_from_member( T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8 )
        : member( x1, x2, x3, x4, x5, x6, x7, x8 ) 
        {}

    template< typename T1, typename T2, typename T3, typename T4, typename T5,
     typename T6, typename T7, typename T8, typename T9 >
    base_from_member( T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
     T9 x9 )
        : member( x1, x2, x3, x4, x5, x6, x7, x8, x9 ) 
        {}

    template< typename T1, typename T2, typename T3, typename T4, typename T5,
     typename T6, typename T7, typename T8, typename T9, typename T10 >
    base_from_member( T1 x1, T2 x2, T3 x3, T4 x4, T5 x5, T6 x6, T7 x7, T8 x8,
     T9 x9, T10 x10 )
        : member( x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 ) 
        {}

};  // boost::base_from_member

}  // namespace boost


#endif  // BOOST_UTILITY_BASE_FROM_MEMBER_HPP
