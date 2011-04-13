// ======================================================================
//
// Copyright (c) 1999,2001 The CGAL Consortium
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
// file          : Threetuple.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL__THREETUPLE_H
#define CGAL__THREETUPLE_H

CGAL_BEGIN_NAMESPACE

template < class T >
class _Threetuple : public Rep
{
public:

  T  e0;
  T  e1;
  T  e2;

  _Threetuple()
  {}

  _Threetuple(const T & a0, const T & a1, const T & a2)
    : e0(a0), e1(a1), e2(a2)
  {}

  ~_Threetuple()
  {}
};

template < class T >
struct Threetuple : public Ref_counted
{
  T  e0;
  T  e1;
  T  e2;

  Threetuple()
  {}

  Threetuple(const T & a0, const T & a1, const T & a2)
    : e0(a0), e1(a1), e2(a2)
  {}
};

template < class T >
class Simple_Threetuple
{
public:
  T  e0;
  T  e1;
  T  e2;

  Simple_Threetuple()
  {}

  Simple_Threetuple(const T & a0, const T & a1, const T & a2)
    : e0(a0), e1(a1), e2(a2)
  {}
};

CGAL_END_NAMESPACE

#endif // CGAL__THREETUPLE_H
