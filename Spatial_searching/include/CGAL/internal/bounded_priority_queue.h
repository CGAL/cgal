// Copyright (c) 2007-10 INRIA (FRANCE).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Gael Guennebaud (gael.guennebaud@inria.fr)




#ifndef CGAL_INTERNAL_BOUNDED_PRIORITY_QUEUE_H
#define CGAL_INTERNAL_BOUNDED_PRIORITY_QUEUE_H

#include <CGAL/license/Spatial_searching.h>


#include <vector>
#include <functional>
#include <algorithm>
#include <boost/next_prior.hpp>

namespace CGAL {
namespace internal{

/**
  * A priority queue with fixed maximum capacity.
  * While the queue has not reached its maximum capacity, elements are
  * inserted as they will be in a heap, the root (top()) being such that
  * Compare(top(),x)=false for any x in the queue.
  * Once the queue is full, trying to insert x in the queue will have no effect if 
  * Compare(x,top())=false. Otherwise, the element at the root of the heap is removed
  * and x is inserted so as to keep the heap property.
  */
template <typename T, typename Compare = std::less<T> >
class bounded_priority_queue
{
public:

  typedef T value_type;
  typedef typename std::vector<value_type>::const_iterator const_iterator;

  bounded_priority_queue(const Compare& comp = Compare())
    : m_comp(comp)
  {}

  bounded_priority_queue(int size, const Compare& comp = Compare())
    : m_count(0), m_data(size), m_comp(comp)
  {}

  /** Sets the max number of elements in the queue */
  void resize(int new_size)
  {
    if (m_data.size()!=new_size)
      m_data.resize(new_size);
    clear();
  }

  /** \returns the number of elements in the queue */
  inline unsigned int size() const { return m_count; }

  /** Removes all elements of the queue. The max size remains unchanged. */
  inline void clear() { m_count = 0; }

  inline bool full() const { return m_count == m_data.size(); }
  inline bool empty() const { return m_count == 0; }

  /** \returns greatest element */
  inline const value_type& top() const { return m_data[0]; }

  inline void insert(const value_type& x)
  {
    value_type* data1 = (&m_data[0]-1);
    if (full())
    {
      if (m_comp(x, top()))
      {
        //insert x in the heap at the correct place,
        //going down in the tree.
        unsigned int j(1), k(2);
        while (k <= m_count)
        {
          value_type* z = &(data1[k]);
          if ((k < m_count) && m_comp(*z, data1[k+1]))
            z = &(data1[++k]);

          if (m_comp(*z, x))
            break;
          data1[j] = *z;
          j = k;
          k = j << 1; //a son of j in the tree
        }
        data1[j] = x;
      }
    }
    else
    {
      //insert element as in a heap
      int i(++m_count), j;
      while (i >= 2)
      {
        j = i >> 1; //father of i in the tree
        value_type& y = data1[j];
        if (m_comp(x, y))
          break;
        data1[i] = y;
        i = j;
      }
      data1[i] = x;
    }
  }

  const_iterator begin() const { return m_data.begin(); }
  const_iterator end() const
  {
    const_iterator res = m_data.begin();
    res += m_count;
    return res;
  }

  void sort()
  {
std::sort(m_data.begin(), boost::next(m_data.begin(),m_count), m_comp);
  }

protected:

  unsigned int m_count;
  std::vector<value_type> m_data;
  Compare m_comp;
};


} } //namespace CGAL::internal

#endif //CGAL_INTERNAL_BOUNDED_PRIORITY_QUEUE_H
