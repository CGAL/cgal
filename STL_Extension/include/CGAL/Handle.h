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
	: PTR(static_cast<Rep*>(0)) {}

    Handle(const Handle& x)
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      PTR = x.PTR;
      PTR->count++;
    }

    ~Handle()
    {
	if ( PTR && (--PTR->count == 0))
	    delete PTR;
    }

    Handle&
    operator=(const Handle& x)
    {
      CGAL_precondition( x.PTR != static_cast<Rep*>(0) );
      x.PTR->count++;
      if ( PTR && (--PTR->count == 0))
	  delete PTR;
      PTR = x.PTR;
      return *this;
    }

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
