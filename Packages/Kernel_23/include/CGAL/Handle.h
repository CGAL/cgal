// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : Handle.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_HANDLE_H
#define CGAL_HANDLE_H

#include <CGAL/Handle_for.h>

CGAL_BEGIN_NAMESPACE

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
    Handle()
	: PTR(static_cast<Rep*>(0)) {}

    Handle(const Handle& x)
    {
      CGAL_kernel_precondition( x.PTR != static_cast<Rep*>(0) );
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
      CGAL_kernel_precondition( x.PTR != static_cast<Rep*>(0) );
      x.PTR->count++;
      if ( PTR && (--PTR->count == 0))
	  delete PTR;
      PTR = x.PTR;
      return *this;
    }

    int
    refs()  const { return PTR->count; }

    friend unsigned long id(const Handle& x);

  protected:
    Rep* PTR;
};

inline
unsigned long
id(const Handle& x)
{ return reinterpret_cast<unsigned long>(x.PTR); }

template < class T >
inline
bool
identical(const T &t1, const T &t2)
{ return id(t1) == id(t2); }

CGAL_END_NAMESPACE

#endif // CGAL_HANDLE_H
