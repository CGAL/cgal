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

class Leda_like_rep
{
    friend class Leda_like_handle;
  protected:
             Leda_like_rep() { count = 1; }
    virtual ~Leda_like_rep() {}

    int      count;
};

class Leda_like_handle
{
  public:
    Leda_like_handle()
	: PTR(static_cast<Leda_like_rep*>(0)) {}

    Leda_like_handle(const Leda_like_handle& x)
    {
      CGAL_kernel_precondition( x.PTR != static_cast<Leda_like_rep*>(0) );
      PTR = x.PTR;
      PTR->count++;
    }

    ~Leda_like_handle()
    {
	if ( PTR && (--PTR->count == 0))
	    delete PTR;
    }

    Leda_like_handle&
    operator=(const Leda_like_handle& x)
    {
      CGAL_kernel_precondition( x.PTR != static_cast<Leda_like_rep*>(0) );
      x.PTR->count++;
      if ( PTR && (--PTR->count == 0))
	  delete PTR;
      PTR = x.PTR;
      return *this;
    }

    int
    refs()  const { return PTR->count; }

    friend unsigned long id(const Leda_like_handle& x);

  protected:
    Leda_like_rep* PTR;
};

inline
unsigned long
id(const Leda_like_handle& x)
{ return reinterpret_cast<unsigned long>(x.PTR); }

template < class T >
inline
bool
identical(const T &t1, const T &t2)
{ return id(t1) == id(t2); }

CGAL_END_NAMESPACE

#if defined(CGAL_USE_LEDA) && !defined(CGAL_NO_LEDA_HANDLE)
#  include <LEDA/basic.h>

CGAL_BEGIN_NAMESPACE

typedef handle_base      Handle;
typedef handle_rep       Rep;

inline
unsigned long
id(const Handle& x)
{ return ID_Number(x); }

CGAL_END_NAMESPACE

#  else

CGAL_BEGIN_NAMESPACE

typedef Leda_like_handle Handle;
typedef Leda_like_rep    Rep;

CGAL_END_NAMESPACE

#endif // defined(CGAL_USE_LEDA) && !defined(CGAL_NO_LEDA_HANDLE)
#endif // CGAL_HANDLE_H
