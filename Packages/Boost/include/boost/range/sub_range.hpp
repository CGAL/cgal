// Boost.Range library
//
//  Copyright Thorsten Ottosen 2003-2004. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
// For more information, see http://www.boost.org/libs/range/
//

#ifndef BOOST_RANGE_SUB_RANGE_HPP
#define BOOST_RANGE_SUB_RANGE_HPP

#include <boost/range/config.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/result_iterator.hpp>
#include <boost/range/size_type.hpp>
#include <boost/range/difference_type.hpp>
#include <boost/assert.hpp>

namespace boost
{
    
    template< class ForwardRange > 
    class sub_range : public iterator_range< BOOST_DEDUCED_TYPENAME range_result_iterator<ForwardRange>::type > 
    {
        typedef BOOST_DEDUCED_TYPENAME range_result_iterator<ForwardRange>::type iterator_t;
        typedef iterator_range< iterator_t  > base;

        typedef BOOST_DEDUCED_TYPENAME base::impl impl;
    public:
        typedef BOOST_DEDUCED_TYPENAME range_value<ForwardRange>::type            value_type;
        typedef BOOST_DEDUCED_TYPENAME range_result_iterator<ForwardRange>::type  iterator;
        typedef BOOST_DEDUCED_TYPENAME range_const_iterator<ForwardRange>::type   const_iterator;
        typedef BOOST_DEDUCED_TYPENAME range_difference<ForwardRange>::type       difference_type;
        typedef BOOST_DEDUCED_TYPENAME range_size<ForwardRange>::type             size_type;

    public:
        sub_range() : base() 
        { }
/*
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))  
        typedef sub_range<ForwardRange> this_type;
        
        sub_range( this_type r ) :
        : base( r )
        { }

        this_type& operator=( this_type r )
        {
            base::operator=( r );
            return *this;
        }
#endif
*/
        template< class ForwardRange2 >
        sub_range( ForwardRange2& r ) : 
            
#if BOOST_WORKAROUND(BOOST_INTEL_CXX_VERSION, <= 800 )
            base( impl::adl_begin( r ), impl::adl_end( r ) )
#else
            base( r )
#endif        
        { }
        
        template< class ForwardRange2 >
        sub_range( const ForwardRange2& r ) : 

#if BOOST_WORKAROUND(BOOST_INTEL_CXX_VERSION, <= 800 )
            base( impl::adl_begin( r ), impl::adl_end( r ) )
#else
            base( r )
#endif                
        { }

        template< class Iter >
        sub_range( Iter first, Iter last ) :
            base( first, last )
        { }
        
        template< class ForwardRange2 >
        sub_range& operator=( ForwardRange2& r )
        {
            base::operator=( r );
            return *this;
        }

        template< class ForwardRange2 >
        sub_range& operator=( const ForwardRange2& r )
        {
            base::operator=( r );
            return *this;
        }
        
    public:
        
        iterator        begin()          { return base::begin(); }
        const_iterator  begin() const    { return base::begin(); }
        iterator        end()            { return base::end();   }
        const_iterator  end() const      { return base::end();   }
        size_type       size() const     { return base::size();  }   

        
    public: // convenience
        value_type& front()
        {
            return base::front();
        }

        const value_type& front() const
        {
            return base::front();
        }

        value_type& back()
        {
            return base::back();
        }

        const value_type& back() const
        {
            return base::back();
        }

        value_type& operator[]( size_type sz )
        {
            return base::operator[](sz);
        }

        const value_type& operator[]( size_type sz ) const
        {
            return base::operator[](sz);
        }

    };

    template< class ForwardRange, class ForwardRange2 >
    inline bool operator==( const sub_range<ForwardRange>& l,
                            const sub_range<ForwardRange2>& r )
    {
        return iterator_range_detail::equal( l, r );
    }

    template< class ForwardRange, class ForwardRange2 >
    inline bool operator!=( const sub_range<ForwardRange>& l,
                            const sub_range<ForwardRange2>& r )
    {
        return !iterator_range_detail::equal( l, r );
    }

    template< class ForwardRange, class ForwardRange2 >
    inline bool operator<( const sub_range<ForwardRange>& l,
                           const sub_range<ForwardRange2>& r )
    {
        return iterator_range_detail::less_than( l, r );
    }


} // namespace 'boost'

#endif
