// Boost.Range library
//
//  Copyright Thorsten Ottosen 2003-2004. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
// For more information, see http://www.boost.org/libs/range/
//

#ifndef BOOST_RANGE_END_HPP
#define BOOST_RANGE_END_HPP

#if defined(_MSC_VER) && (_MSC_VER >= 1200)
# pragma once
#endif

#include <boost/range/config.hpp>

#ifdef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
#include <boost/range/detail/end.hpp>
#else

#include <boost/range/detail/implementation_help.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/const_iterator.hpp>

namespace boost 
{ 
namespace range_detail
{

        //////////////////////////////////////////////////////////////////////
        // primary template
        //////////////////////////////////////////////////////////////////////
        
        template< typename C >
        inline BOOST_DEDUCED_TYPENAME range_const_iterator<C>::type
        end( const C& c )
        {
            return c.end();
        }
                
        template< typename C >
        inline BOOST_DEDUCED_TYPENAME range_iterator<C>::type
        end( C& c )
        {
            return c.end();
        }
      
        //////////////////////////////////////////////////////////////////////
        // pair
        //////////////////////////////////////////////////////////////////////

        template< typename Iterator >
        inline Iterator end( const std::pair<Iterator,Iterator>& p )
        {
            return p.second;
        }
        
        template< typename Iterator >
        inline Iterator end( std::pair<Iterator,Iterator>& p )
        {
            return p.second;
        }
        
        //////////////////////////////////////////////////////////////////////
        // array
        //////////////////////////////////////////////////////////////////////

        template< typename T, std::size_t sz >
        inline const T* end( const T (&array)[sz] )
        {
            return range_detail::array_end<T,sz>( array ); 
        }
        
        template< typename T, std::size_t sz >
        inline T* end( T (&array)[sz] )
        {
            return range_detail::array_end<T,sz>( array ); 
        }

        //////////////////////////////////////////////////////////////////////
        // string
        //////////////////////////////////////////////////////////////////////

#if BOOST_WORKAROUND(__MWERKS__, <= 0x3204 ) || BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
// CW up to 9.3 and borland have troubles with function ordering
        inline char* end( char* s )
        {
            return range_detail::str_end( s );
        }

        inline wchar_t* end( wchar_t* s )
        {
            return range_detail::str_end( s );
        }

        inline const char* end( const char* s )
        {
            return range_detail::str_end( s );
        }

        inline const wchar_t* end( const wchar_t* s )
        {
            return range_detail::str_end( s );
        }
#else
        inline char* end( char*& s )
        {
            return range_detail::str_end( s );
        }

        inline wchar_t* end( wchar_t*& s )
        {
            return range_detail::str_end( s );
        }

        inline const char* end( const char*& s )
        {
            return range_detail::str_end( s );
        }

        inline const wchar_t* end( const wchar_t*& s )
        {
            return range_detail::str_end( s );
        }
#endif
        
} // namespace 'range_detail'

template< class T >
inline BOOST_DEDUCED_TYPENAME range_iterator<T>::type end( T& r )
{
    return range_detail::end( r );
}

template< class T >
inline BOOST_DEDUCED_TYPENAME range_const_iterator<T>::type end( const T& r )
{
    return range_detail::end( r );
}



#if BOOST_WORKAROUND(__MWERKS__, <= 3003 ) || BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
// BCB and CW are not able to overload pointer when class overloads are also available.
template<>
inline range_const_iterator<const char*>::type end<const char*>( const char*& r )
{
    return range_detail::str_end( r );
}

template<>
inline range_const_iterator<const wchar_t*>::type end<const wchar_t*>( const wchar_t*& r )
{
    return range_detail::str_end( r );
}

#endif

} // namespace 'boost'



#endif // BOOST_NO_FUNCTION_TEMPLATE_ORDERING


namespace boost
{
    template< class T >
    inline BOOST_DEDUCED_TYPENAME range_const_iterator<T>::type
    const_end( const T& r )
    {
        return end( r );
    }
}

#endif
