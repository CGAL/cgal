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
#include <atomic>
#include <CGAL/Handle_for.h>
#include <CGAL/assertions.h>

namespace CGAL {

class Rep
{
    friend class Handle;
  protected:
             Rep() : count(1) { }
    virtual ~Rep() {}

    std::atomic_int count;
};

class Handle
{
  public:

    typedef std::ptrdiff_t Id_type ;

    Handle() noexcept
        : PTR(static_cast<Rep*>(0)) {}

    // FIXME: if the precondition throws in a noexcept function, the program terminates
    Handle(const Handle& x) noexcept
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      PTR = x.PTR;
      //CGAL_assume (PTR->count > 0);
      incref();
    }

    Handle(Handle&& x) noexcept : PTR(x.PTR) { x.PTR = 0; }

    ~Handle() { reset(); }

    Handle&
    operator=(const Handle& x) noexcept
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      x.incref();
      reset();
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

    void reset()
    {
      if (PTR)
      {
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
        PTR=0;
      }
    }

    void incref() const noexcept
    {
      PTR->count.fetch_add(1, std::memory_order_relaxed);
    }

    int
    refs()  const noexcept { return PTR->count; }

    Id_type id() const noexcept { return PTR - static_cast<Rep*>(0); }

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
