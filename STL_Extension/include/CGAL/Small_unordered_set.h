// Copyright (c) 2020 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_SMALL_UNORDERED_SET_H
#define CGAL_SMALL_UNORDERED_SET_H

#include <array>
#include <unordered_set>

namespace CGAL
{

/*
  This is a very rudimentary structure. It is far from being a full
  "small unordered set", but it is a starting point.

  For the moment, its only feature is the insertion of elements +
  check that they were indeed inserted, while avoiding instantiating a
  `std::unordered_set` as long as the number of elements remain small.

  In practice:

  - if the number of elements in the set is lower than MaxSize,
    elements are just appended to a `std::array<Key, MaxSize>` and the
    unicity test is done element by element, in linear time

  - when the number of elements exceed MaxSize, a
    `std::unordered_set<Key>` is instanciated, all the elements of the
    array are inserted in it and from that point the container behaves
    like a `std::unordered_set`

  For this structure to be a true "small unordered set", special
  iterators should be created to make the switch from array to
  unordered set completely transparent to the user, and all other
  usual member functions should be introduced. So far, this is not
  needed and thus not done.
*/
template<typename Key, std::size_t MaxSize>
class Small_unordered_set
{
  using Array = std::array<Key, MaxSize>;
  using Set   = std::unordered_set<Key>;

  Array m_array;
  std::unique_ptr<Set> m_set;
  std::size_t m_size = 0;

public:

  Small_unordered_set() { }

  Small_unordered_set (const Small_unordered_set& other)
    : m_size (other.m_size)
  {
    if (other.m_set)
      m_set = std::make_unique<Set>(*other.m_set);
    else
      m_array = other.m_array;
  }

  Small_unordered_set (Small_unordered_set&& other)
    : m_size (other.m_size)
  {
    if (other.m_set)
      m_set = std::move(other.m_set);
    else
      m_array = std::move(other.m_array);
  }

  Small_unordered_set& operator= (const Small_unordered_set& other)
  {
    m_size = other.m_size;
    if (other.m_set)
      m_set = std::make_unique<Set>(*other.m_set);
    else
      m_array = other.m_array;
  }

  Small_unordered_set& operator= (Small_unordered_set&& other)
  {
    m_size = other.m_size;
    if (other.m_set)
      m_set = std::move(other.m_set);
    else
      m_array = std::move(other.m_array);
  }

  bool insert (const Key& key)
  {
    if (m_size != MaxSize)
    {
      for (std::size_t i = 0; i < m_size; ++ i)
        if (m_array[i] == key)
          return false;
      m_array[m_size ++] = key;
      return true;
    }

    if (!m_set)
    {
      m_set = std::make_unique<Set>();
      m_set->reserve (MaxSize + 1);
      for (const Key& a : m_array)
        m_set->insert(a);
    }

    return m_set->insert(key).second;
  }

};

}


#endif // CGAL_SMALL_UNORDERED_SET_H
