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
// release_date  : 2000, December 10
// 
// source        : Tuple.fw
// file          : Twotuple.h
// package       : Kernel_basic (3.17)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.17
// revision_date : 10 Dec 2000 
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

CGAL_END_NAMESPACE

#endif // CGAL__TWOTUPLE_H
