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
    
    Handle()
      : PTR{static_cast<Rep*>(0)} {}

    Handle(const Handle& x)
    {
      CGAL_precondition( x.PTR.p != static_cast<Rep*>(0) );
      PTR.p = x.PTR.p;
      PTR.p->count++;
    }

    ~Handle()
    {
	if ( PTR.p && (--PTR.p->count == 0))
	    delete PTR.p;
    }

    Handle&
    operator=(const Handle& x)
    {
      CGAL_precondition( x.PTR.p != static_cast<Rep*>(0) );
      x.PTR.p->count++;
      if ( PTR.p && (--PTR.p->count == 0))
	  delete PTR.p;
      PTR.p = x.PTR.p;
      return *this;
    }

    void reset()
    {
      if (PTR.p)
      {
        if (--PTR.p->count==0)
          delete PTR.p;
        PTR.p=0;
      }
    }

    int
    refs()  const { return PTR.p->count; }

    Id_type id() const { return PTR.p - static_cast<Rep*>(0); }

    bool identical(const Handle& h) const { return PTR.p == h.PTR.p; }

    void*  for_compact_container() const { return PTR.vp; }
    void*& for_compact_container() { return PTR.vp; }

  protected:

  union {
    Rep* p;
    void* vp;
  } PTR;
};

//inline Handle::Id_type id(const Handle& x) { return x.id() ; }

inline bool identical(const Handle &h1, const Handle &h2) { return h1.identical(h2); }

} //namespace CGAL

#endif // CGAL_HANDLE_H
