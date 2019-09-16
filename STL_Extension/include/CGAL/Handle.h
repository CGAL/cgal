// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_HANDLE_H
#define CGAL_HANDLE_H

#include <cstddef>
#include <CGAL/Handle_for.h>
#include <CGAL/assertions.h>

namespace CGAL {

class Rep
{
    friend class Handle;
  protected:
    Rep() : count(1) {}
    virtual ~Rep() {}

#if defined(CGAL_HANDLE_FOR_USE_ATOMIC) && ! defined(CGAL_NO_ATOMIC)
    CGAL::cpp11::atomic<unsigned int> count;
#elif CGAL_HANDLE_FOR_USE_BOOST_ATOMIC_COUNTER
    boost::detail::atomic_count count;
#else // no atomic
    unsigned int count;
#endif // no atomic

};

class Handle
{
  public:
    
    typedef std::ptrdiff_t Id_type ;
    
    Handle()
	: PTR(static_cast<Rep*>(0)) {}

    Handle(const Handle& x)
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      CGAL_assume (x.PTR->count > 0);
      PTR = x.PTR;
#if defined(CGAL_HANDLE_FOR_USE_ATOMIC) && ! defined(CGAL_NO_ATOMIC)
      PTR->count.fetch_add(1, CGAL::cpp11::memory_order_relaxed);
#else // not CGAL::cpp11::atomic
      ++(PTR->count);
#endif // not CGAL::cpp11::atomic
    }

    ~Handle()
    {
      if ( PTR ) {
        CGAL_assume (PTR->count > 0);
#if defined(CGAL_HANDLE_FOR_USE_ATOMIC) && ! defined(CGAL_NO_ATOMIC)
        if (PTR->count.fetch_sub(1, CGAL::cpp11::memory_order_release) == 1) {
          CGAL::cpp11::atomic_thread_fence(CGAL::cpp11::memory_order_acquire);
          delete PTR;
        }
#else // not CGAL::cpp11::atomic
        if (--(PTR->count) == 0) {
          delete PTR;
        }
#endif // not CGAL::cpp11::atomic
      }
    }

    Handle&
    operator=(const Handle& x)
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      CGAL_assume (x.PTR->count > 0);
#if defined(CGAL_HANDLE_FOR_USE_ATOMIC) && ! defined(CGAL_NO_ATOMIC)
      x.PTR->count.fetch_add(1, CGAL::cpp11::memory_order_relaxed);
#else // not CGAL::cpp11::atomic
      ++(x.PTR->count);
#endif // not CGAL::cpp11::atomic
      if ( PTR ) {
        CGAL_assume (PTR->count > 0);
#if defined(CGAL_HANDLE_FOR_USE_ATOMIC) && ! defined(CGAL_NO_ATOMIC)
        if (PTR->count.fetch_sub(1, CGAL::cpp11::memory_order_release) == 1) {
          CGAL::cpp11::atomic_thread_fence(CGAL::cpp11::memory_order_acquire);
          delete PTR;
        }
#else // not CGAL::cpp11::atomic
        if (--(PTR->count) == 0) {
          delete PTR;
        }
#endif // not CGAL::cpp11::atomic
      }
      PTR = x.PTR;
      return *this;
    }

    int
    refs()  const { return PTR->count; }

    Id_type id() const { return PTR - static_cast<Rep*>(0); }

    bool identical(const Handle& h) const { return PTR == h.PTR; }

  protected:
    Rep* PTR;
};

//inline Handle::Id_type id(const Handle& x) { return x.id() ; }

inline bool identical(const Handle &h1, const Handle &h2) { return h1.identical(h2); }

} //namespace CGAL

#endif // CGAL_HANDLE_H
