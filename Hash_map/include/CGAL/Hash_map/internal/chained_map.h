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
#ifndef CGAL_HASH_MAP_INTERNAL_CHAINED_MAP_H
#define CGAL_HASH_MAP_INTERNAL_CHAINED_MAP_H

#include <CGAL/memory.h>
#include <CGAL/assertions.h>
#include <iostream>
#include <limits>
#include <type_traits>
#include <utility>

namespace CGAL {

namespace internal {

template <typename T, typename Allocator = CGAL_ALLOCATOR(T) > class chained_map;
template <typename T> class chained_map_elem;

template <typename T>
class chained_map_elem
{
  template<typename T2, typename Alloc> friend class chained_map;
  std::size_t k; T i;
  chained_map_elem<T>*  succ;
};

template <typename T, typename Allocator>
class chained_map
{
   static constexpr std::size_t nullkey = (std::numeric_limits<std::size_t>::max)();

   chained_map_elem<T>* table;
   chained_map_elem<T>* table_end;
   chained_map_elem<T>* freelist;
   std::size_t table_size;
   std::size_t table_size_1;

   typedef std::allocator_traits<Allocator> Allocator_traits;
   typedef typename Allocator_traits::template rebind_alloc<chained_map_elem<T> > allocator_type;

   allocator_type alloc;
   std::size_t reserved_size;
   T def;

public:
   T& xdef() { return def; }
   const T& cxdef() const { return def; }
private:
   void init_inf(T& x)   const { x = def; }


   chained_map_elem<T>*  HASH(std::size_t x)  const
   { return table + (x & table_size_1);  }

   void init_table(std::size_t n);
   void rehash();

   inline void insert(std::size_t x, T y);

   void destroy(chained_map_elem<T>* item)
   {
     typedef std::allocator_traits<allocator_type> Allocator_type_traits;
     Allocator_type_traits::destroy(alloc,item);
   }

public:
   static constexpr std::size_t min_size = 32;
   static constexpr std::size_t default_size = 512;
   typedef chained_map_elem<T>*  Item;

   std::size_t index(Item it) const { return it->k; }
   T&            inf(Item it) const { return it->i; }

   chained_map(std::size_t n = default_size, const T& d = T());
   chained_map(const chained_map<T, Allocator>& D);
   chained_map& operator=(const chained_map<T, Allocator>& D);
   chained_map(chained_map<T, Allocator>&& D)
     noexcept(std::is_nothrow_move_constructible_v<Allocator> && std::is_nothrow_move_constructible_v<T>);
   chained_map& operator=(chained_map<T, Allocator>&& D)
     noexcept(std::is_nothrow_move_assignable_v<Allocator> && std::is_nothrow_move_assignable_v<T>);

   void reserve(std::size_t n);
   void clear();
   ~chained_map()
   {
     if(!table)
       return;
     for (Item item = table ; item != table_end ; ++item)
       destroy(item);
     alloc.deallocate(table, table_end - table);
   }

   T& access(Item p, std::size_t x);
   T& access(std::size_t x);
   Item lookup(std::size_t x) const;
   void statistics() const;
};

template <typename T, typename Allocator>
inline T& chained_map<T, Allocator>::access(std::size_t x)
{
  if(!table)
    init_table(reserved_size);

  Item p = HASH(x);

  if ( p->k == x ) {
     return p->i;
  }
  else {
    if ( p->k == nullkey ) {
      p->k = x;
      init_inf(p->i);  // initializes p->i to xdef
      return p->i;
    } else
      return access(p,x);
  }
}

template <typename T, typename Allocator>
void chained_map<T, Allocator>::init_table(std::size_t n)
{
  std::size_t t = min_size;
  while (t < n) t <<= 1;

  table_size = t;
  table_size_1 = t-1;
  table = alloc.allocate(t + t/2);
  for (std::size_t i = 0 ; i < t + t/2 ; ++i){
    std::allocator_traits<allocator_type>::construct(alloc,table + i);
  }

  freelist = table + t;
  table_end = table + t + t/2;

  for (Item p = table; p < freelist; ++p)
  { p->succ = nullptr;
    p->k = nullkey;
  }
}


template <typename T, typename Allocator>
inline void chained_map<T, Allocator>::insert(std::size_t x, T y)
{ Item q = HASH(x);
  if ( q->k == nullkey ) {
    q->k = x;
    q->i = y;
  } else {
    freelist->k = x;
    freelist->i = y;
    freelist->succ = q->succ;
    q->succ = freelist++;
  }
}


template <typename T, typename Allocator>
void chained_map<T, Allocator>::rehash()
{
  chained_map_elem<T>* old_table = table;
  chained_map_elem<T>* old_table_end = table_end;

  Item old_table_mid = table + table_size;

  init_table(2*table_size);

  Item p;

  for(p = old_table; p < old_table_mid; ++p)
  { std::size_t x = p->k;
    if ( x != nullkey ) // list p is non-empty
    { Item q = HASH(x);
      q->k = x;
      q->i = p->i;
    }
  }

  while (p < old_table_end)
  { std::size_t x = p->k;
    insert(x,p->i);
    ++p;
  }

  for (Item item = old_table ; item != old_table_end ; ++item)
    destroy(item);
  alloc.deallocate(old_table, old_table_end - old_table);
}


template <typename T, typename Allocator>
T& chained_map<T, Allocator>::access(Item p, std::size_t x)
{
  Item q = p->succ;
  while (q && q->k != x) q = q->succ;
  if (q)
  {
    return q->i;
  }

  // index x not present, insert it

  if (freelist == table_end)   // table full: rehash
  { rehash();
    p = HASH(x);
  }

  if (p->k == nullkey)
  { p->k = x;
    init_inf(p->i);  // initializes p->i to xdef
    return p->i;
  }

  q = freelist++;
  q->k = x;
  init_inf(q->i);    // initializes q->i to xdef
  q->succ = p->succ;
  p->succ = q;
  return q->i;
}


template <typename T, typename Allocator>
chained_map<T, Allocator>::chained_map(std::size_t n, const T& d)
  : table(nullptr), reserved_size(n), def(d)
{
}


template <typename T, typename Allocator>
chained_map<T, Allocator>::chained_map(const chained_map<T, Allocator>& D)
{
  init_table(D.table_size);

  for(Item p = D.table; p < D.freelist; ++p)
  { if (p->k != nullkey || p >= D.table + D.table_size)
    { insert(p->k,p->i);
      //D.copy_inf(p->i);  // see chapter Implementation
    }
  }
}
template <typename T, typename Allocator>
chained_map<T, Allocator>::chained_map(chained_map<T, Allocator>&& D)
  noexcept(std::is_nothrow_move_constructible_v<Allocator> && std::is_nothrow_move_constructible_v<T>)
  : table(std::exchange(D.table, nullptr))
  , table_end(std::exchange(D.table_end, nullptr))
  , freelist(std::exchange(D.freelist, nullptr))
  , table_size(std::exchange(D.table_size, 0))
  , table_size_1(std::exchange(D.table_size_1, 0))
  , alloc(std::move(D.alloc))
  , reserved_size(std::exchange(D.reserved_size, 0))
  , def(std::move(D.def))
{}

template <typename T, typename Allocator>
chained_map<T, Allocator>& chained_map<T, Allocator>::operator=(const chained_map<T, Allocator>& D)
{
  clear();

  init_table(D.table_size);

  for(Item p = D.table; p < D.freelist; ++p)
  { if (p->k != nullkey || p >= D.table + D.table_size)
    { insert(p->k,p->i);
      //copy_inf(p->i);    // see chapter Implementation
    }
  }
  return *this;
}

template <typename T, typename Allocator>
chained_map<T, Allocator>& chained_map<T, Allocator>::operator=(chained_map<T, Allocator>&& D)
  noexcept(std::is_nothrow_move_assignable_v<Allocator> && std::is_nothrow_move_assignable_v<T>)
{
  clear();

  table = std::exchange(D.table, nullptr);
  table_end = std::exchange(D.table_end, nullptr);
  freelist = std::exchange(D.freelist, nullptr);
  table_size = std::exchange(D.table_size, 0);
  table_size_1 = std::exchange(D.table_size_1, 0);
  alloc = std::move(D.alloc);
  reserved_size = std::exchange(D.reserved_size, 0);
  def = std::move(D.def);

  return *this;
}

template <typename T, typename Allocator>
void chained_map<T, Allocator>::reserve(std::size_t n)
{
  CGAL_assertion(!table);
  reserved_size = n;
}

template <typename T, typename Allocator>
void chained_map<T, Allocator>::clear()
{
  if(!table)
    return;

  for (Item item = table ; item != table_end ; ++item)
    destroy(item);
  alloc.deallocate(table, table_end - table);

  table = nullptr;
}

template <typename T, typename Allocator>
typename chained_map<T, Allocator>::Item
chained_map<T, Allocator>::lookup(std::size_t x) const
{
  if(!table)
    return nullptr;

  Item p = HASH(x);
  while (p && p->k != x)
  { p = p->succ; }
  return p;
}

template <typename T, typename Allocator>
void chained_map<T, Allocator>::statistics() const
{ std::cout << "table_size: " << table_size <<"\n";
  std::size_t n = 0;
  for (Item p = table; p < table + table_size; ++p)
     if (p ->k != nullkey) ++n;
  std::size_t used_in_overflow = freelist - (table + table_size );
  n += used_in_overflow;
  std::cout << "number of entries: " << n << "\n";
  std::cout << "fraction of entries in first position: " <<
               ((double) (n - used_in_overflow))/n <<"\n";
  std::cout << "fraction of empty lists: " <<
               ((double) (n - used_in_overflow))/table_size<<"\n";
}

} // namespace internal
} //namespace CGAL

#endif // CGAL_HASH_MAP_INTERNAL_CHAINED_MAP_H
