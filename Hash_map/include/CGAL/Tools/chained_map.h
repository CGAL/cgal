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
   static constexpr std::size_t nullptrKEY = (std::numeric_limits<std::size_t>::max)();

   chained_map_elem<T>* table;
   chained_map_elem<T>* table_end;
   chained_map_elem<T>* free;
   std::size_t table_size;
   std::size_t table_size_1;
   T def;

   typedef std::allocator_traits<Allocator> Allocator_traits;
   typedef typename Allocator_traits::template rebind_alloc<chained_map_elem<T> > allocator_type;

   allocator_type alloc;
   std::size_t reserved_size;

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
   typedef chained_map_elem<T>*  chained_map_item;
   typedef chained_map_item item;

   std::size_t index(chained_map_item it) const { return it->k; }
   T&            inf(chained_map_item it) const { return it->i; }

   chained_map(std::size_t n = default_size, const T& d = T());
   chained_map(const chained_map<T, Allocator>& D);
   chained_map& operator=(const chained_map<T, Allocator>& D);

   void reserve(std::size_t n);
   void clear();
   ~chained_map()
   {
     if(!table)
       return;
     for (chained_map_item item = table ; item != table_end ; ++item)
       destroy(item);
     alloc.deallocate(table, table_end - table);
   }

   T& access(chained_map_item p, std::size_t x);
   T& access(std::size_t x);
   chained_map_item lookup(std::size_t x) const;
   void statistics() const;
};

template <typename T, typename Allocator>
inline T& chained_map<T, Allocator>::access(std::size_t x)
{
  if(!table)
    init_table(reserved_size);

  chained_map_item p = HASH(x);

  if ( p->k == x ) {
     return p->i;
  }
  else {
    if ( p->k == nullptrKEY ) {
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

  free = table + t;
  table_end = table + t + t/2;

  for (chained_map_item p = table; p < free; p++)
  { p->succ = nullptr;
    p->k = nullptrKEY;
  }
}


template <typename T, typename Allocator>
inline void chained_map<T, Allocator>::insert(std::size_t x, T y)
{ chained_map_item q = HASH(x);
  if ( q->k == nullptrKEY ) {
    q->k = x;
    q->i = y;
  } else {
    free->k = x;
    free->i = y;
    free->succ = q->succ;
    q->succ = free++;
  }
}


template <typename T, typename Allocator>
void chained_map<T, Allocator>::rehash()
{
  chained_map_elem<T>* old_table = table;
  chained_map_elem<T>* old_table_end = table_end;
  chained_map_elem<T>* old_free = free;

  chained_map_item old_table_mid = table + table_size;

  init_table(2*table_size);

  chained_map_item p;

  for(p = old_table; p < old_table_mid; p++)
  { std::size_t x = p->k;
    if ( x != nullptrKEY ) // list p is non-empty
    { chained_map_item q = HASH(x);
      q->k = x;
      q->i = p->i;
    }
  }

  while (p < old_table_end)
  { std::size_t x = p->k;
    insert(x,p->i);
    p++;
  }

  for (chained_map_item item = old_table ; item != old_table_end ; ++item)
    destroy(item);
  alloc.deallocate(old_table, old_table_end - old_table);
}


template <typename T, typename Allocator>
T& chained_map<T, Allocator>::access(chained_map_item p, std::size_t x)
{
  chained_map_item q = p->succ;
  while (q && q->k != x) q = q->succ;
  if (q)
  {
    return q->i;
  }

  // index x not present, insert it

  if (free == table_end)   // table full: rehash
  { rehash();
    p = HASH(x);
  }

  if (p->k == nullptrKEY)
  { p->k = x;
    init_inf(p->i);  // initializes p->i to xdef
    return p->i;
  }

  q = free++;
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

  for(chained_map_item p = D.table; p < D.free; p++)
  { if (p->k != nullptrKEY || p >= D.table + D.table_size)
    { insert(p->k,p->i);
      //D.copy_inf(p->i);  // see chapter Implementation
    }
  }
}

template <typename T, typename Allocator>
chained_map<T, Allocator>& chained_map<T, Allocator>::operator=(const chained_map<T, Allocator>& D)
{
  clear();

  init_table(D.table_size);

  for(chained_map_item p = D.table; p < D.free; p++)
  { if (p->k != nullptrKEY || p >= D.table + D.table_size)
    { insert(p->k,p->i);
      //copy_inf(p->i);    // see chapter Implementation
    }
  }
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

  for (chained_map_item item = table ; item != table_end ; ++item)
    destroy(item);
  alloc.deallocate(table, table_end - table);

  table = nullptr;
}

template <typename T, typename Allocator>
typename chained_map<T, Allocator>::chained_map_item
chained_map<T, Allocator>::lookup(std::size_t x) const
{
  if(!table)
    return nullptr;

  chained_map_item p = HASH(x);
  while (p && p->k != x)
  { p = p->succ; }
  return p;
}

template <typename T, typename Allocator>
void chained_map<T, Allocator>::statistics() const
{ std::cout << "table_size: " << table_size <<"\n";
  std::size_t n = 0;
  for (chained_map_item p = table; p < table + table_size; p++)
     if (p ->k != nullptrKEY) n++;
  std::size_t used_in_overflow = free - (table + table_size );
  n += used_in_overflow;
  std::cout << "number of entries: " << n << "\n";
  std::cout << "fraction of entries in first position: " <<
               ((double) (n - used_in_overflow))/n <<"\n";
  std::cout << "fraction of empty lists: " <<
               ((double) (n - used_in_overflow))/table_size<<"\n";
}

} // namespace internal
} //namespace CGAL

#endif // CGAL_CHAINED_MAP_H
