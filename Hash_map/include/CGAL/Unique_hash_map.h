// Copyright (c) 1997-2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
//                 Lutz Kettner <kettner@inf.ethz.ch>

#ifndef CGAL_UNIQUE_HASH_MAP_H
#define CGAL_UNIQUE_HASH_MAP_H

#include <CGAL/disable_warnings.h>

#include <CGAL/config.h>
#include <CGAL/memory.h>
#include <CGAL/Handle_hash_function.h>
#include <CGAL/Hash_map/internal/chained_map.h>
#include <cstddef>

namespace CGAL {

template <class Key_, class Data_,
          class UniqueHashFunction = Handle_hash_function,
          class Allocator_ = CGAL_ALLOCATOR(Data_) >
class Unique_hash_map {
public:
    typedef Key_                                     Key;
    typedef Data_                                    Data;
    typedef UniqueHashFunction                       Hash_function;
    typedef Allocator_                               Allocator;

    // STL compliance
    typedef Key_                                     key_type;
    typedef Data_                                    data_type;
    typedef UniqueHashFunction                       hasher;

    typedef Unique_hash_map<Key,Data,Hash_function,Allocator> Self;

private:
    typedef internal::chained_map<Data, Allocator>   Map;
    typedef typename Map::Item                       Item;

private:
    Hash_function  m_hash_function;
    Map            m_map;

    template <class It, class Iterator_category>
    void reserve_impl(It, It, Iterator_category)
    {}

    template <class It>
    void reserve_impl(It b, It e, std::forward_iterator_tag)
    {
      m_map.reserve(std::distance(b,e));
    }

public:

    Unique_hash_map() = default;

    Unique_hash_map( const Data& deflt, std::size_t table_size = Map::default_size)
        : m_map(table_size, deflt)
    {}

    Unique_hash_map( const Data& deflt,
                     std::size_t table_size,
                     const Hash_function& fct)
        : m_hash_function(fct), m_map( table_size, deflt)
    {}

    Unique_hash_map( Key first1, Key beyond1, Data first2) {
        insert( first1, beyond1, first2);
    }
    Unique_hash_map( Key first1, Key beyond1, Data first2,
                     const Data& deflt,
                     std::size_t table_size   = 1,
                     const Hash_function& fct = Hash_function())
    : m_hash_function(fct), m_map(table_size, deflt) {
        insert( first1, beyond1, first2);
    }

    void reserve(std::size_t n)
    { m_map.reserve(n); }

    Data default_value() const { return m_map.cxdef(); }

    Hash_function  hash_function() const { return m_hash_function; }

    void clear() { m_map.clear(); }

    void clear( const Data& deflt) {
        m_map.clear();
        m_map.xdef() = deflt; }

    bool is_defined( const Key& key) const {
        return m_map.lookup( m_hash_function(key)) != 0;
    }

    const Data& operator[]( const Key& key) const {
        Item p = m_map.lookup( m_hash_function(key));
        if ( p != 0 )
            return m_map.inf(p);
        return m_map.cxdef();
    }

    Data& operator[]( const Key& key) {
        return m_map.access( m_hash_function(key));
    }

    Data insert( Key first1, Key beyond1, Data first2) {
        reserve_impl(first1, beyond1, typename std::iterator_traits<Key>::iterator_category());
        for ( ; first1 != beyond1; (++first1, ++first2)) {
            operator[]( first1) = first2;
        }
        return first2;
    }

    void statistics() const { m_map.statistics(); }
};



} //namespace CGAL

namespace boost {
  template <typename UniquePairAssociativeContainer>
  class associative_property_map;

  struct lvalue_property_map_tag;

  template <typename KeyType, typename ValueType,
            typename HashFunction, typename Allocator>
  class associative_property_map<CGAL::Unique_hash_map<KeyType, ValueType,
                                                       HashFunction, Allocator> >
  {
    typedef CGAL::Unique_hash_map<KeyType, ValueType, HashFunction, Allocator> C;

  public:
    typedef KeyType key_type;
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef lvalue_property_map_tag category;

    associative_property_map() : m_c(0) { }
    associative_property_map(C& c) : m_c(&c) { }

    reference operator[](const key_type& k) const {
      return (*m_c)[k];
    }

    friend
    reference
    get(const associative_property_map<C>& uhm, const key_type& key)
    {
      return uhm[key];
    }

    friend
    void
    put(associative_property_map<C>& uhm, const key_type& key, const value_type& val)
    {
      uhm[key] = val;
    }

  private:
    C* m_c;
  };

  template <typename KeyType, typename ValueType,
            typename HashFunction, typename Allocator>
  associative_property_map<CGAL::Unique_hash_map<KeyType, ValueType,
                                                 HashFunction, Allocator> >
  make_assoc_property_map(CGAL::Unique_hash_map<KeyType, ValueType,
                                                HashFunction, Allocator>& c)
  {
    return associative_property_map<CGAL::Unique_hash_map<KeyType, ValueType,
                                                          HashFunction, Allocator> >(c);
  }

}

#include <CGAL/enable_warnings.h>

#endif // CGAL_UNIQUE_HASH_MAP_H
// EOF
