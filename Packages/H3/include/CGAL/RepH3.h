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
// file          : RepH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_REPH3_H
#define CGAL_REPH3_H

CGAL_BEGIN_NAMESPACE

template <class NT>
class RepH3 : public Ref_counted
{
public:
  NT  e0;
  NT  e1;
  NT  e2;
  NT  e3;

  RepH3()
   // : e0(NT(42)), e1(NT(42)), e2(NT(42)), e3(NT(1))
  {}
  RepH3(const NT& a0, const NT& a1, const NT& a2, const NT& a3)
    : e0(a0), e1(a1), e2(a2), e3(a3)
  {}

  NT    hx() { return e0; }
  NT    hy() { return e1; }
  NT    hz() { return e2; }
  NT    hw() { return e3; } // homogenizing component
};


CGAL_END_NAMESPACE


#endif // CGAL_REPH3_H
