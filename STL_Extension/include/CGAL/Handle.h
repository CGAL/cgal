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
#include <CGAL/Handle_for.h>
#include <CGAL/assertions.h>

namespace CGAL {

class Rep
{
    friend class Handle;
  protected:
             Rep() { count = 1; }
    virtual ~Rep() {}

    int      count;
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
      CGAL_assume (PTR->count > 0);
      PTR->count++;
    }

    ~Handle()
    {
        if ( PTR && (--PTR->count == 0))
            delete PTR;
    }

    Handle&
    operator=(const Handle& x) noexcept
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      x.PTR->count++;
      if ( PTR && (--PTR->count == 0))
          delete PTR;
      PTR = x.PTR;
      return *this;
    }

    friend void swap(Handle& a, Handle& b) noexcept { std::swap(a.PTR, b.PTR); }

    void reset()
    {
      if (PTR)
      {
        if (--PTR->count==0)
          delete PTR;
        PTR=0;
      }
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
