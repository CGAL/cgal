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
#include <CGAL/unordered_flat_map.h>
#include <cstddef>
#include <iostream>
#include <iterator>

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
    typedef std::size_t internal_key;

    struct identity_hash {
        std::size_t operator()(internal_key x) const noexcept { return x; }
    };

    typedef std::pair<const internal_key, Data> entry;
    typedef std::allocator_traits<Allocator> traits;
    typedef typename traits::template rebind_alloc<entry> entry_allocator;

    typedef CGAL::unordered_flat_map<internal_key, Data, identity_hash,
        std::equal_to<internal_key>, entry_allocator>  Map;

    Hash_function  m_hash_function;
    Map            m_map;
    Data           m_default_value;

    template <class It, class Iterator_category>
    void reserve_impl(It, It, Iterator_category)
    {}

    template <class It>
    void reserve_impl(It b, It e, std::forward_iterator_tag)
    {
      m_map.reserve(std::distance(b,e));
    }

public:
    static constexpr std::size_t default_size = 8;

    Unique_hash_map(const Data& deflt = Data(),
                    std::size_t table_size = default_size,
                    const Hash_function& fct = Hash_function())
        : m_hash_function(fct), m_map(table_size), m_default_value(deflt) {
    }

    Unique_hash_map(Key first1, Key beyond1, Data first2,
                    const Data& deflt = Data(),
                    std::size_t table_size = default_size,
                    const Hash_function& fct = Hash_function())
    : m_hash_function(fct), m_map(table_size), m_default_value(deflt) {
        insert( first1, beyond1, first2);
    }

    void reserve(std::size_t n) {
        m_map.reserve(n);
    }

    Data default_value() const {
        return m_default_value;
    }

    Hash_function hash_function() const {
        return m_hash_function;
    }

    void clear() {
        m_map.clear();
    }

    void clear(const Data& deflt) {
        m_map.clear();
        m_default_value = deflt;
    }

    bool is_defined(const Key& key) const {
        return m_map.find(m_hash_function(key)) != m_map.end();
    }

    bool contains(const Key& key) const {
        return m_map.contains(m_hash_function(key));
    }

    std::size_t erase(const Key& key) {
        return m_map.erase(m_hash_function(key));
    }

    const Data& operator[](const Key& key) const {
        auto it = m_map.find(m_hash_function(key));
        if (it != m_map.end()) return it->second;
        return m_default_value;
    }

    Data& operator[](const Key& key) {
        return m_map.try_emplace(
            m_hash_function(key), m_default_value).first->second;
    }

    Data insert( Key first1, Key beyond1, Data first2) {
        reserve_impl(first1, beyond1, typename std::iterator_traits<Key>::iterator_category());
        for ( ; first1 != beyond1; (++first1, ++first2)) {
            operator[]( first1) = first2;
        }
        return first2;
    }

    void statistics() const {
        std::cout << "CGAL::unordered_flat_map: size=" << m_map.size() << "\n";
    }
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
