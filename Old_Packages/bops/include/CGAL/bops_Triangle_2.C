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
// file          : include/CGAL/bops_Triangle_2.C
// package       : bops (2.2)
// source        : include/CGAL/bops_Triangle_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL_BOPS_TRIANGLE_2_C
#define CGAL_BOPS_TRIANGLE_2_C

//#define CGAL__BOPS_USE_STANDARD_CODE
   // if CGAL__BOPS_USE_STANDARD_CODE is defined, then the
   // code for simple polygons will be used for computing 
   // BOPS on iso-rectangles and triangles.

#include <CGAL/bops_simple_polygons_2.h>
#include <CGAL/bops_Container_Polygon_2.h>

#include<CGAL/bops_Triangle_2.h>
#include <CGAL/bops_assertions.h>

CGAL_BEGIN_NAMESPACE 

/*
#include <CGAL/Triangle_2_Triangle_2_intersection.h>
bool do_intersect( const Triangle_2<R>& A,
                        const Triangle_2<R>& B);
Object intersection( const Triangle_2<R>& A,
                        const Triangle_2<R>& B)
see 'CGAL/Triangle_2_Triangle_2_intersection.h'
*/
 

template < class R, class OutputIterator >
OutputIterator Union( const Triangle_2<R>& A,
                           const Triangle_2<R>& B,
		           OutputIterator result)
{
  CGAL_bops_precondition_msg( !A.is_degenerate(),
                              "Triangle_2<R> A is degenerated");
  CGAL_bops_precondition_msg( !B.is_degenerate(),
                              "Triangle_2<R> B is degenerated");

  Bops_default_I<R> default_traits;
  typedef typename Bops_default_I<R>::Input_polygon_container input;
  input cA, cB;
  back_insert_iterator<input> inA(cA), inB(cB);
  for(int i= 0; i<3; i++) { *inA++ = A[i]; *inB++ = B[i]; }

  return Union( cA.begin(), cA.end(), cB.begin(), cB.end(),
	 default_traits, result );
}
 
 
template < class R, class OutputIterator >
OutputIterator difference( const Triangle_2<R>& A,
		                const Triangle_2<R>& B,
		                OutputIterator result)
{
  CGAL_bops_precondition_msg( !A.is_degenerate(),
                              "Triangle_2<R> A is degenerated");
  CGAL_bops_precondition_msg( !B.is_degenerate(),
                              "Triangle_2<R> B is degenerated");

  Bops_default_I<R> default_traits;
  typedef typename Bops_default_I<R>::Input_polygon_container input;
  input cA, cB;
  back_insert_iterator<input> inA(cA), inB(cB);
  for(int i= 0; i<3; i++) { *inA++= A[i]; *inB++= B[i]; }

  return difference( cA.begin(), cA.end(), cB.begin(), cB.end(),
	 default_traits, result );
} 
 
CGAL_END_NAMESPACE

#endif // CGAL_BOPS_TRIANGLE_2_C
