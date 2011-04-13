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
// file          : Twotuple.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL__TWOTUPLE_H
#define CGAL__TWOTUPLE_H

CGAL_BEGIN_NAMESPACE

template < class T >
class _Twotuple : public Rep
{
public:
  T  e0;
  T  e1;

  _Twotuple()
  {}

  _Twotuple(const T & a0, const T &a1)
  : e0(a0), e1(a1)
  {}

  ~_Twotuple()
  {}
};

template < class T >
class Twotuple : public Ref_counted
{
public:
  T  e0;
  T  e1;

  Twotuple()
  {}

  Twotuple(const T & a0, const T &a1) : e0(a0), e1(a1)
  {}
};

template < class T >
class Simple_Twotuple
{
public:
  T  e0, e1;

  Simple_Twotuple()
  {}

  Simple_Twotuple(const T & a0, const T &a1)
  : e0(a0), e1(a1)
  {}
};

CGAL_END_NAMESPACE

#endif // CGAL__TWOTUPLE_H
