// ======================================================================
//
// Copyright (c) 2000,2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Kernel_d/Line_d.C
// package       : Kernel_d
// maintainer    : Michael Seel <Michael.Seel@mpi-sb.mpg.de>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Seel
// coordinator   : MPI Saarbruecken (Susan.Hert@mpi-sb.mpg.de)
// ======================================================================

#ifndef CGAL_LINE_D_C
#define CGAL_LINE_D_C
CGAL_BEGIN_NAMESPACE

template <class R> 
Line_d<R> Segment_d<R>::supporting_line() const
{ CGAL_assertion_msg((!is_degenerate()), 
  "Segment_d::supporting_line(): degenerate segment cannot be converted.");
  return Line_d<R>(Base(*this)); 
} 

template <class R>
Line_d<R> Ray_d<R>::supporting_line() const
{ return Line_d<R>(Base(*this)); } 

CGAL_END_NAMESPACE
#endif //CGAL_LINE_D_C

