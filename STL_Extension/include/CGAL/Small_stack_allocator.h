// Copyright (c) 2020  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_SMALL_STACK_ALLOCATOR_H
#define CGAL_SMALL_STACK_ALLOCATOR_H

#include <CGAL/Profile_counter.h>

#include <cstdlib>


namespace CGAL
{

namespace internal
{

template <typename T, std::size_t MaxSize>
class Small_stack_allocator_pool
{
  T m_pool[MaxSize];
  std::size_t m_next;

public:

  Small_stack_allocator_pool ()
    : m_next(0) { }
  ~Small_stack_allocator_pool() { }

  // Forbid copy/assignment
  Small_stack_allocator_pool (const Small_stack_allocator_pool&) = delete;
  Small_stack_allocator_pool& operator= (const Small_stack_allocator_pool&) = delete;

  T* allocate (std::size_t n)
  {
    CGAL_BRANCH_PROFILER("stack allocations / total allocations of Small_stack_allocator", prof);
    if (m_next + n < MaxSize)
    {
      CGAL_BRANCH_PROFILER_BRANCH(prof);
      T* out = m_pool + m_next;
      m_next += n;
      return out;
    }

    return new T[n];
  }

  void deallocate (T* p, std::size_t n)
  {
    if (m_pool <= p && p < m_pool + MaxSize)
    {
      std::size_t pos = static_cast<std::size_t>(p - m_pool);
      if (pos + n == m_next)
        m_next = pos;
    }
    else
      delete[] p;
  }

};

} // namespace internal

template <typename T, std::size_t MaxSize>
class Small_stack_allocator
{
public:

  using value_type = T;

  template <typename T2, std::size_t MaxSize2>
  friend class Small_stack_allocator;

  template <typename T2>
  struct rebind
  {
    using other = Small_stack_allocator<T2, MaxSize>;
  };

private:

  using Pool = internal::Small_stack_allocator_pool<T, MaxSize>;
  Pool m_pool;

public:

  Small_stack_allocator() { }

  // Forbid copy/assignment
  Small_stack_allocator (const Small_stack_allocator&) = delete;
  Small_stack_allocator& operator=  (const Small_stack_allocator&) = delete;

  T* allocate (std::size_t n)
  {
    return m_pool.allocate(n);
  }

  void deallocate (T* p, std::size_t n)
  {
    m_pool.deallocate (p, n);
  }
};


} // namespace CGAL



#endif // CGAL_SMALL_STACK_ALLOCATOR_H
