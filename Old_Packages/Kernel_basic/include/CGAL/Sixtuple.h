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
// file          : Sixtuple.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL__SIXTUPLE_H
#define CGAL__SIXTUPLE_H

CGAL_BEGIN_NAMESPACE

template < class T >
class _Sixtuple : public Rep
{
public:

  T  e0;
  T  e1;
  T  e2;
  T  e3;
  T  e4;
  T  e5;

  _Sixtuple()
   {
   }
  _Sixtuple(const T & a0, const T & a1, const T & a2,
                 const T & a3, const T & a4, const T & a5)
    : e0(a0), e1(a1), e2(a2), e3(a3), e4(a4), e5(a5)
  {}

  ~_Sixtuple()
  {}
};

template < class T >
class Sixtuple : public Ref_counted
{
public:

  T  e0;
  T  e1;
  T  e2;
  T  e3;
  T  e4;
  T  e5;

  Sixtuple()
  {}

  Sixtuple(const T & a0, const T & a1, const T & a2,
           const T & a3, const T & a4, const T & a5)
    : e0(a0), e1(a1), e2(a2), e3(a3), e4(a4), e5(a5)
  {}
};

template < class T >
class Simple_Sixtuple
{
public:

  T  e0;
  T  e1;
  T  e2;
  T  e3;
  T  e4;
  T  e5;

  Simple_Sixtuple()
  {}

  Simple_Sixtuple(const T & a0, const T & a1, const T & a2,
                  const T & a3, const T & a4, const T & a5)
    : e0(a0), e1(a1), e2(a2), e3(a3), e4(a4), e5(a5)
  {}
};

CGAL_END_NAMESPACE

#endif // CGAL__SIXTUPLE_H
