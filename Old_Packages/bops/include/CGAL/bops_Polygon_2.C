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
// file          : include/CGAL/bops_Polygon_2.C
// package       : bops (2.2)
// source        : include/CGAL/bops_Polygon_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_BOPS_POLYGON_2_C
#define CGAL_BOPS_POLYGON_2_C

#include <CGAL/bops_Container_Polygon_2.h>
#include <CGAL/bops_Polygon_2.h>
#include <CGAL/bops_assertions.h>
#include <CGAL/bops_Convex_Polygon_2.h>

CGAL_BEGIN_NAMESPACE

template < class R, class Container >
bool do_intersect(
     const Polygon_2<Polygon_traits_2<R>, Container>& A,
     const Polygon_2<Polygon_traits_2<R>, Container>& B)
{
  CGAL_bops_precondition_msg(A.is_simple(),
                             "Polygon_2<R> A is not simple");
  CGAL_bops_precondition_msg(B.is_simple(),
                             "Polygon_2<R> B is not simple");
  CGAL_bops_precondition_msg(A.is_counterclockwise_oriented(),
                             "Polygon_2<R> A is not counterclockwise oriented");
  CGAL_bops_precondition_msg(B.is_counterclockwise_oriented(),
                             "Polygon_2<R> B is not counterclockwise oriented");

  Bops_default_I<R> default_traits;
  return do_intersect(
         A.vertices_begin(), A.vertices_end(),
         B.vertices_begin(), B.vertices_end(),
	 default_traits );
}

 
template < class R, class Container, class OutputIterator >
OutputIterator intersection(
       const Polygon_2<Polygon_traits_2<R>, Container>& A,
       const Polygon_2<Polygon_traits_2<R>, Container>& B,
       OutputIterator result)
{
  CGAL_bops_precondition_msg(A.is_simple(),
                             "Polygon_2<R> A is not simple");
  CGAL_bops_precondition_msg(B.is_simple(),
                             "Polygon_2<R> B is not simple");
  CGAL_bops_precondition_msg(A.is_counterclockwise_oriented(),
                             "Polygon_2<R> A is not counterclockwise oriented");
  CGAL_bops_precondition_msg(B.is_counterclockwise_oriented(),
                             "Polygon_2<R> B is not counterclockwise oriented");


  if( A.is_convex() && B.is_convex() ) {
     Polygon_2<Polygon_traits_2<R>, Container> C, AA(A), BB(B);
     C= Convex_Intersection(AA, BB);
     result++ = make_object(C);
     return result;
  }

  Bops_default_I<R> default_traits;
  return intersection(
         A.vertices_begin(), A.vertices_end(),
         B.vertices_begin(), B.vertices_end(),
	 default_traits, result );
}


template < class R, class Container, class OutputIterator >
OutputIterator Union(
      const Polygon_2<Polygon_traits_2<R>, Container>& A,
      const Polygon_2<Polygon_traits_2<R>, Container>& B,
      OutputIterator result)
{
  CGAL_bops_precondition_msg(A.is_simple(),
                             "Polygon_2<R> A is not simple");
  CGAL_bops_precondition_msg(B.is_simple(),
                             "Polygon_2<R> B is not simple");
  CGAL_bops_precondition_msg(A.is_counterclockwise_oriented(),
                             "Polygon_2<R> A is not counterclockwise oriented");
  CGAL_bops_precondition_msg(B.is_counterclockwise_oriented(),
                             "Polygon_2<R> B is not counterclockwise oriented");

  Bops_default_I<R> default_traits;
  return Union(
         A.vertices_begin(), A.vertices_end(),
         B.vertices_begin(), B.vertices_end(),
	 default_traits, result );
}
 
template < class R, class Container, class OutputIterator >
OutputIterator difference(
	       const Polygon_2<Polygon_traits_2<R>, Container>& A,
	       const Polygon_2<Polygon_traits_2<R>, Container>& B,
	       OutputIterator result)
{
  CGAL_bops_precondition_msg(A.is_simple(),
                             "Polygon_2<R> A is not simple");
  CGAL_bops_precondition_msg(B.is_simple(),
                             "Polygon_2<R> B is not simple");
  CGAL_bops_precondition_msg(A.is_counterclockwise_oriented(),
                             "Polygon_2<R> A is not counterclockwise oriented");
  CGAL_bops_precondition_msg(B.is_counterclockwise_oriented(),
                             "Polygon_2<R> B is not counterclockwise oriented");

  Bops_default_I<R> default_traits;
  return difference(
         A.vertices_begin(), A.vertices_end(),
         B.vertices_begin(), B.vertices_end(),
	 default_traits, result );
} 

CGAL_END_NAMESPACE

#endif // CGAL_BOPS_POLYGON_2_C
