// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 01
//
// file          : include/CGAL/bops_Iso_rectangle_2.h
// package       : bops (2.2)
// source        : include/CGAL/bops_Iso_rectangle_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_BOPS_ISO_RECTANGLE_2_H
#define CGAL_BOPS_ISO_RECTANGLE_2_H

#include <CGAL/bops_traits_2.h>
#include <CGAL/bops_assertions.h>


#include <CGAL/Iso_rectangle_2_Iso_rectangle_2_intersection.h>

CGAL_BEGIN_NAMESPACE

/*
bool 
do_intersect( const Iso_rectangle_2<R>& A,
                   const Iso_rectangle_2<R>& B);
*/

template < class R, class OutputIterator >
inline OutputIterator
intersection( const Iso_rectangle_2<R>& A,
                   const Iso_rectangle_2<R>& B,
                   OutputIterator object_it) {
  CGAL_bops_precondition_msg( !A.is_degenerate(),
                              "Iso_rectangle_2<R> A is degenerated");
  CGAL_bops_precondition_msg( !B.is_degenerate(),
                              "Iso_rectangle_2<R> B is degenerated");

  *object_it++= intersection(A,B);
  return object_it;
}


template < class R, class OutputIterator >
OutputIterator
Union(        const Iso_rectangle_2<R>& A,
                   const Iso_rectangle_2<R>& B,
                   OutputIterator object_it);
 
 
template < class R, class OutputIterator >
OutputIterator
difference(   const Iso_rectangle_2<R>& A,
                   const Iso_rectangle_2<R>& B,
                   OutputIterator object_it);
 
CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/bops_Iso_rectangle_2.C>
#endif

#endif // CGAL_BOPS_ISO_RECTANGLE_2_H
