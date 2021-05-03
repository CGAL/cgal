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
// Author(s)     : Courtesy of LEDA
#ifndef CGAL_CHAINED_MAP_H
#define CGAL_CHAINED_MAP_H

#include <CGAL/memory.h>
#include <iostream>

namespace CGAL {

namespace internal {

template <typename T, typename Allocator = CGAL_ALLOCATOR(T)>
class chained_map
{
  struct element_type {
    std::size_t key;
    T value;
    element_type* next = nullptr;
  };

public:
  typedef element_type* item;
  typedef std::size_t key_type;
  typedef std::size_t size_type;

  static constexpr size_type min_size = 32;

  chained_map() : reserved_size(min_size)
  {
  }

  chained_map(size_type n) : reserved_size(n)
  {
  }

  ~chained_map()
  {
#ifdef CHAINED_MAP_DEBUG
    size_type n = number_of_entries();
    if (reserved_size < n)
      std::cout << "reserved: " << reserved_size << " entries: " << n << std::endl;
#endif
    destroy(table);
  }

  void reserve(size_type n)
  {
    reserved_size = n;
  }

  T& xdef()
  {
    return stop.value;
  }

  const T& cxdef() const
  {
    return stop.value;
  }

  T& inf(item it) const
  {
    return it->value;
  }

  void clear()
  {
    if (!table.begin)
      return;

    for (item it = table.begin + 1; it < table.free; ++it)
      if (it->key != nullptr_key || it >= table.begin + table.size)
        it->value = T();

    destroy(table);
    table.begin = nullptr;
  }

  T& access(key_type key)
  {
    if (!table.begin)
      init_table();

    item it = hash(key);

    if (it->key == key) {
      old_key = key;
      return it->value;
    }

    if (it->key == nullptr_key) {
      it->key = key;
      it->value = stop.value; // initializes p->value to xdef
      old_key = key;
      return it->value;
    }

    return access(it, key);
  }

  item lookup(key_type key) const
  {
    if (!table.begin)
      return nullptr;

    item it = hash(key);
    const_cast<key_type&>(stop.key) = key; // cast away const
    while (it->key != key) {
      it = it->next;
    }
    return (it == &stop) ? nullptr : it;
  }

  size_type number_of_entries()
  {
    size_type n = 0;
    if (!table.begin)
      return n;

    for (item it = table.begin + 1; it < table.begin + table.size; ++it)
      if (it->key != nullptr_key)
        ++n;

    return n;
  }

  void statistics() const
  {
    std::cout << "table_size: " << table.size << std::endl;
    size_type n = number_of_entries();
    size_type used_in_overflow = table.free - (table.begin + table.size);
    n += used_in_overflow;
    std::cout << "number of entries: " << n << std::endl;
    std::cout << "fraction of entries in first position: "
              << ((double)(n - used_in_overflow)) / n << std::endl;
    std::cout << "fraction of empty lists: "
              << ((double)(n - used_in_overflow)) / table.size << std::endl;
  }

private:
  struct table_type {
    item begin = nullptr;
    item end = nullptr;
    item free = nullptr;
    size_type size;
    size_type size_1;
  };

  typedef std::allocator_traits<Allocator> allocator_traits;
  typedef typename allocator_traits::template rebind_alloc<element_type> allocator_type;
  typedef std::allocator_traits<allocator_type> allocator_type_traits;

  item hash(key_type key) const
  {
    return table.begin + (key & table.size_1);
  }

  void init_table()
  {
    init_table(reserved_size);
  }

  void init_table(size_type n)
  {
    size_type t = min_size;
    while (t < n)
      t <<= 1;

    table.size = t;
    table.size_1 = t - 1;
    size_type s = t + t / 2;
    table.begin = alloc.allocate(s);
    for (size_type i = 0; i < s; ++i)
      allocator_type_traits::construct(alloc, table.begin + i);

    table.free = table.begin + t;
    table.end = table.begin + s;

    for (item it = table.begin; it < table.free; ++it) {
      it->next = &stop;
      it->key = nullptr_key;
    }

    table.begin->key = non_nullptr_key;
  }

  void rehash()
  {
    table_type old_table = table;

    item old_table_mid = table.begin + table.size;

    init_table(2 * table.size);

    item it;
    for (it = old_table.begin + 1; it < old_table_mid; ++it) {
      key_type key = it->key;
      if (key != nullptr_key) // list p is non-empty
      {
        item q = hash(key);
        q->key = key;
        q->value = it->value;
      }
    }

    for (; it < old_table.end; ++it) {
      key_type x = it->key;
      insert(x, it->value);
    }

    // delete old table
    table_type cur_table = table;
    table = old_table;
    old_table.begin = nullptr;
    T p = access(old_key);
    destroy(table);
    table = cur_table;
    access(old_key) = p;
  }

  void insert(key_type key, T value)
  {
    item q = hash(key);
    if (q->key == nullptr_key) {
      q->key = key;
      q->value = value;
    } else {
      item f = table.free;
      f->key = key;
      f->value = value;
      f->next = q->next;
      q->next = table.free++;
    }
  }

  T& access(item it, key_type key)
  {
    stop.key = key;
    item q = it->next;
    while (q->key != key)
      q = q->next;

    if (q != &stop) {
      old_key = key;
      return q->value;
    }

    // index x not present, insert it

    if (table.free == table.end) { // table full: rehash
      rehash();
      it = hash(key);
    }

    if (it->key == nullptr_key) {
      it->key = key;
      it->value = stop.value; // initializes p->value to xdef
      return it->value;
    }

    q = table.free++;
    q->key = key;
    q->value = stop.value; // initializes q->value to xdef
    q->next = it->next;
    it->next = q;
    return q->value;
  }

  void destroy(const table_type& t)
  {
    if (!t.begin)
      return;

    for (item it = t.begin; it != t.end; ++it)
      destroy(it);

    alloc.deallocate(t.begin, t.end - t.begin);
  }

  void destroy(item it)
  {
    allocator_type_traits::destroy(alloc, it);
  }

  static constexpr key_type nullptr_key = 0;
  static constexpr key_type non_nullptr_key = 1;

  element_type stop;
  table_type table;
  size_type reserved_size;
  key_type old_key;
  allocator_type alloc;
};

} // namespace internal
} // namespace CGAL

#endif // CGAL_CHAINED_MAP_H
