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
// file          : Fourtuple.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL__FOURTUPLE_H
#define CGAL__FOURTUPLE_H

CGAL_BEGIN_NAMESPACE

template < class T >
class _Fourtuple : public Rep
{
public:
  T  e0;
  T  e1;
  T  e2;
  T  e3;

  _Fourtuple()
  {}
  _Fourtuple(const T & a0, const T & a1, const T & a2, const T & a3)
    : e0(a0), e1(a1), e2(a2), e3(a3)
  {}

  ~_Fourtuple()
  {}
};

template < class T >
class Fourtuple : public Ref_counted
{
public:
  T  e0;
  T  e1;
  T  e2;
  T  e3;

  Fourtuple()
  {}

  Fourtuple(const T & a0, const T & a1, const T & a2, const T & a3)
    : e0(a0), e1(a1), e2(a2), e3(a3)
  {}
};

template < class T >
class Simple_Fourtuple
{
public:
  T  e0;
  T  e1;
  T  e2;
  T  e3;

  Simple_Fourtuple()
  {}

  Simple_Fourtuple(const T & a0, const T & a1, const T & a2, const T & a3)
    : e0(a0), e1(a1), e2(a2), e3(a3)
  {}
};

CGAL_END_NAMESPACE

#endif // CGAL__FOURTUPLE_H
