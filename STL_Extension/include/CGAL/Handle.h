// Copyright (c) 1999
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_HANDLE_H
#define CGAL_HANDLE_H

#include <cstddef>
#include <cstdint>
#include <atomic>
#include <CGAL/Handle_for.h>
#include <CGAL/assertions.h>

namespace CGAL {

class Rep
{
    friend class Handle;
  protected:
    Rep(int count = 1)
      : count(count)
    {}
    virtual ~Rep() {}

    std::atomic_int count;
};


class Handle
{
  public:

    typedef std::ptrdiff_t Id_type ;

    Handle() noexcept
        : PTR(static_cast<Rep*>(0)) {}

    Handle(const Handle& x) noexcept(!(CGAL_PRECONDITIONS_ENABLED || CGAL_ASSERTIONS_ENABLED))
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      PTR = x.PTR;
      //CGAL_assume (PTR->count > 0);
      incref();
    }

    Handle(Handle&& x) noexcept : PTR(x.PTR) { x.PTR = 0; }

    ~Handle() { reset(); }

    Handle&
    operator=(const Handle& x) noexcept(!CGAL_PRECONDITIONS_ENABLED)
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      x.incref();
      if(PTR) decref(); // not reset() in case this==&x
      PTR = x.PTR;
      return *this;
    }

    Handle&
    operator=(Handle&& x) noexcept
    {
      swap(*this, x);
      return *this;
    }

    friend void swap(Handle& a, Handle& b) noexcept { std::swap(a.PTR, b.PTR); }

  private:
    void incref() const noexcept
    {
      if (is_currently_single_threaded()) {
        PTR->count.store(PTR->count.load(std::memory_order_relaxed) + 1, std::memory_order_relaxed);
      } else {
        PTR->count.fetch_add(1, std::memory_order_relaxed);
      }
    }

    void decref()
    {
      if (is_currently_single_threaded()) {
        auto c = PTR->count.load(std::memory_order_relaxed);
        if (c == 1)
          delete PTR;
        else
          PTR->count.store(c - 1, std::memory_order_relaxed);
      } else {
      // TSAN does not support fences :-(
#if !defined __SANITIZE_THREAD__ && !__has_feature(thread_sanitizer)
        if (PTR->count.load(std::memory_order_relaxed) == 1
            || PTR->count.fetch_sub(1, std::memory_order_release) == 1) {
          std::atomic_thread_fence(std::memory_order_acquire);
#else
        if (PTR->count.fetch_sub(1, std::memory_order_acq_rel) == 1) {
#endif
          delete PTR;
        }
      }
    }

  public:
    void reset()
    {
      if (PTR)
      {
        decref();
        PTR=0;
      }
    }

    int
    refs()  const noexcept { return PTR->count.load(std::memory_order_relaxed); }

    Id_type id() const noexcept { return static_cast<Id_type>(reinterpret_cast<std::intptr_t>(static_cast<void*>(PTR)) / sizeof(Rep)); }

    bool identical(const Handle& h) const noexcept { return PTR == h.PTR; }

    void * for_compact_container() const { return PTR; }
    void for_compact_container(void* p) { PTR = static_cast<Rep*>(p); }
  protected:
    Rep* PTR;
};

//inline Handle::Id_type id(const Handle& x) { return x.id() ; }

inline bool identical(const Handle &h1, const Handle &h2) noexcept { return h1.identical(h2); }

} //namespace CGAL

#endif // CGAL_HANDLE_H
