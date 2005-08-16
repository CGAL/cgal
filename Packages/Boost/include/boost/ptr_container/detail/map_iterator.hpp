//
// Boost.Pointer Container
//
//  Copyright Thorsten Ottosen 2003-2005. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
// For more information, see http://www.boost.org/libs/ptr_container/
//

#ifndef BOOST_PTR_CONTAINER_MAP_ITERATOR_HPP
#define BOOST_PTR_CONTAINER_MAP_ITERATOR_HPP

#if defined(_MSC_VER) && (_MSC_VER >= 1200)
# pragma once
#endif

#include <boost/config.hpp>
#include <cassert>
#include <utility>

namespace boost
{
    //namespace ptr_container_detail
    //{
        template
        < 
            typename I, // original iterator 
            typename K, // key type
            typename V  // return value type of operator*()
        > 
        class ptr_map_iterator
        {
            I iter_;
            typedef K              key_type;
            
        public:
            typedef std::ptrdiff_t difference_type;
            typedef V              value_type;
            typedef V*             pointer;
            typedef V&             reference;
            typedef                std::bidirectional_iterator_tag  iterator_category;        
            
        public:
            ptr_map_iterator()                                  {}
            ptr_map_iterator( const I& i ) : iter_( i )         {}
            
            template< class MutableIterator, class K2, class V2 >
            ptr_map_iterator( const ptr_map_iterator<MutableIterator,K2,V2>& r ) 
             : iter_(r.base())
            { }
            
            V& operator*() const
            {
                return *static_cast<V*>( iter_->second );
            }

            V* operator->() const
            {
                return static_cast<V*>( iter_->second );
            }
            
            ptr_map_iterator& operator++()
            {
                ++iter_;
                return *this;
            }

            ptr_map_iterator operator++(int)
            {
                ptr_map_iterator res = *this;
                ++iter_;
                return res;
            }

            ptr_map_iterator& operator--()
            {
                --iter_;
                return *this;
            }

            ptr_map_iterator operator--(int)
            {
                ptr_map_iterator res = *this;
                --iter_;
                return res;

            }

            I base() const
            {
                return iter_;
            }

            key_type key() const
            {
                return iter_->first;
            }

       }; // class 'ptr_map_iterator'


       
       template< class I, class K, class V, class I2, class K2, class V2 >
       inline bool operator==( const ptr_map_iterator<I,K,V>& l, 
                               const ptr_map_iterator<I2,K2,V2>& r )
       {
           return l.base() == r.base();
       }


       
       template< class I, class K, class V, class I2, class K2, class V2 >
       inline bool operator!=( const ptr_map_iterator<I,K,V>& l, 
                               const ptr_map_iterator<I2,K2,V2>& r )
       {
           return l.base() != r.base();
       }
}

#endif
