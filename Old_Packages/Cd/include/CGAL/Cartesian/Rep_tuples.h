// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/Reptuples.
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_CARTESIAN_REP_TUPLES_H
#define CGAL_CARTESIAN_REP_TUPLES_H

// Old tuples deriving from Rep.
// Only used by the old d-dim cartesian kernel.

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

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_REP_TUPLES_H
