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
// file          : include/CGAL/bops_Triangle_2.h
// package       : bops (2.2)
// source        : include/CGAL/bops_Triangle_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_BOPS_TRIANGLE2_H
#define CGAL_BOPS_TRIANGLE2_H


#include <CGAL/bops_traits_2.h>
#include <CGAL/Triangle_2.h>

CGAL_BEGIN_NAMESPACE

/* NOT IMPLEMENTED YET
#include <CGAL/Triangle_2_Triangle_2_intersection.h>
bool 
do_intersect( const Triangle_2<R>& A,
                   const Triangle_2<R>& B);
OutputIterator
intersection( const Triangle_2<R>& A,
                   const Triangle_2<R>& B,
                   OutputIterator object_it);
{
  CGAL_bops_precondition_msg( !A.is_degenerate(),
                              "Triangle_2<R> A is degenerated");
  CGAL_bops_precondition_msg( !B.is_degenerate(),
                              "Triangle_2<R> B is degenerated");

  *object_it++= intersection(A,B);
  return object_it;
}
*/


template < class R, class OutputIterator >
OutputIterator
Union(        const Triangle_2<R>& A,
                   const Triangle_2<R>& B,
                   OutputIterator object_it);
 
 
 
template < class R, class OutputIterator >
OutputIterator
difference(   const Triangle_2<R>& A,
                   const Triangle_2<R>& B,
                   OutputIterator list_of_objects_it);
 
CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/bops_Triangle_2.C>
#endif

#endif // CGAL_BOPS_TRIANGLE2_H
