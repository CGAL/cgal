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
// file          : solve.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_SOLVE_H
#define CGAL_SOLVE_H

CGAL_BEGIN_NAMESPACE


template <class FT>
void solve (const FT &a1, const FT &a2, const FT &a3,
                 const FT &b1, const FT &b2, const FT &b3,
                 const FT &c1, const FT &c2, const FT &c3,
                 const FT &d1, const FT &d2, const FT &d3,
                 FT &x, FT &y, FT &z)
{
  FT denom = b2*c1*a3-b1*c2*a3+c3*b1*a2+b3*c2*a1-c1*b3*a2-b2*c3*a1;

  x = - (b2*c3*d1-b2*c1*d3+c1*b3*d2+b1*c2*d3-c3*b1*d2-b3*c2*d1)/denom;

  z = (b2*d1*a3-b2*a1*d3+b1*a2*d3-b1*d2*a3-d1*b3*a2+a1*b3*d2)/denom;

  y = (a2*c3*d1-a2*c1*d3-c2*d1*a3+c2*a1*d3+d2*c1*a3-d2*c3*a1)/denom;
}


// this is for a parabola c1, c2, c3 are equal to 1
template <class FT>
void solve_quadratic (const FT &a1, const FT &a2, const FT &a3,
                           const FT &b1, const FT &b2, const FT &b3,
                           const FT &d1, const FT &d2, const FT &d3,
                           FT &x, FT &y, FT &z)
{
  FT denom = b2*a3-b1*a3+b1*a2+b3*a1-b3*a2-b2*a1;

  x = - (b2*d1-b2*d3+b3*d2+b1*d3-b1*d2-b3*d1)/denom;

  z = (b2*d1*a3-b2*a1*d3+b1*a2*d3-b1*d2*a3-d1*b3*a2+a1*b3*d2)/denom;

  y = (a2*d1-a2*d3-d1*a3+a1*d3+d2*a3-d2*a1)/denom;
}


CGAL_END_NAMESPACE

#endif // CGAL_SOLVE_H
