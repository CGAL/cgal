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

    // If chunk to allocate does not exceed remaining stack size, use stack
    if (m_next + n < MaxSize)
    {
      CGAL_BRANCH_PROFILER_BRANCH(prof);
      T* out = m_pool + m_next;
      m_next += n;
      return out;
    }

    // Else, fallback to `new` allocation
    return static_cast<T*>(::operator new(n * sizeof(T)));
  }

  void deallocate (T* p, std::size_t n)
  {
    // If pointer is part of the pool, nothing to do
    if (m_pool <= p && p < m_pool + MaxSize)
    {
      std::size_t pos = static_cast<std::size_t>(p - m_pool);

      // If pointer was the last one allocated, we can use again this
      // chunk of memory
      if (pos + n == m_next)
        m_next = pos;
    }
    // Else, delete it
    else
      ::operator delete(p);
  }

};

} // namespace internal

/*
  The small stack allocator holds a statically allocated array of type
  T. When allocation is required, it first returns pointers from that
  statically allocated space. When allocation exceeds this small size,
  it goes to the "normal" mode and allocates with `new`.

  When deallocating, the stack can be used again if the last element
  allocated on it is deallocated (in the absence of the free list, it
  is difficult to do better than that). Keep in mind that this
  allocator is mainly designed for containers that you know are going
  to do few allocations/deallocations.
*/
template <typename T, std::size_t MaxSize>
class Small_stack_allocator
{
public:

  using pointer =               T*;
  using const_pointer =   const T*;
  using reference =             T&;
  using const_reference = const T&;
  using value_type =            T ;

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

  // Do not copy pool
  Small_stack_allocator(const Small_stack_allocator&) { }
  Small_stack_allocator& operator= (const Small_stack_allocator&)
  {
    return Small_stack_allocator();
  }

  T* allocate (std::size_t n)
  {
    return m_pool.allocate(n);
  }

  void deallocate (T* p, std::size_t n)
  {
    m_pool.deallocate (p, n);
  }

  // Functions construct and destroy are not used by modern compilers
  // and are deprecated from the Allocator concept in C++17. They are
  // included to comply with older compilers and should be removed at
  // some point.

  template <typename U, typename ... Args>
  void construct (U* u, Args&& ... args)
  {
    ::new (static_cast<void*>(u))
      U (std::forward<Args> (args)...);
  }

  template <typename U>
  void destroy (U* u)
  {
    u->~U();
  }

};

} // namespace CGAL

#endif // CGAL_SMALL_STACK_ALLOCATOR_H
