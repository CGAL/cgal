// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : 
// source        : chull_traits.lw
// revision      : 2.3  
// revision_date : 01 Feb 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================

#ifndef _TEST_CLS_CHULL_TRAITS_3_C
#define _TEST_CLS_CHULL_TRAITS_3_C

#include <CGAL/basic.h>
#include <CGAL/dd_geo/chull_traits_3.h>
#include <CGAL/dd_geo/chull.h>
#include <vector>

template <class R>
bool
__test_cls_chull_traits_3( const R& )
{
  typedef R                                    RepCls;
  typedef typename R::RT                       RT;
  typedef CGAL::chull_traits_3< RepCls >       ChullTraits;
  typedef typename ChullTraits::POINT          Point;
#ifdef CGAL_CFG_NO_DEFAULT_TEMPLATE_ARGUMENTS
  typedef typename ChullTraits::PLANE          Plane;
  typedef chull< ChullTraits, Point, Plane>    ChullType;
#else
  typedef chull< ChullTraits >                 ChullType;
#endif // CGAL_CFG_NO_DEFAULT_TEMPLATE_ARGUMENTS

  std::vector<Point> V( 8, Point() );
  V[0] = Point( RT(12), RT( 8), RT( 4), RT(2) );
  V[1] = Point( RT(45), RT(-5), RT(20), RT(10) );
  V[2] = Point( RT(-3), RT(19), RT(18), RT(5) );
  V[3] = Point( RT(16), RT(12), RT(-8), RT(4) );
  V[4] = Point( RT(12), RT(27), RT(-23), RT(3) );
  V[5] = Point( RT( 7), RT(46), RT(-12), RT(-4) );
  V[6] = Point( RT(40), RT(55), RT(-25), RT(5) );
  V[7] = Point( RT(-2), RT(-9), RT( 81), RT(6) );

  ChullType CH(3); 
  typename std::vector<Point>::iterator  it;
  for (it = V.begin(); it != V.end(); ++it)  CH.insert(*it);

  V[0] = Point( RT(12), RT( 8), RT( 0), RT(2) );
  V[1] = Point( RT(45), RT(-5), RT( 0), RT(10) );
  V[2] = Point( RT(-3), RT(19), RT( 0), RT(5) );
  V[3] = Point( RT(16), RT(12), RT( 0), RT(4) );
  V[4] = Point( RT(12), RT(27), RT( 0), RT(3) );
  V[5] = Point( RT( 7), RT(46), RT( 0), RT(-4) );
  V[6] = Point( RT(40), RT(55), RT( 0), RT(5) );
  V[7] = Point( RT(-2), RT(-9), RT( 0), RT(6) );

  ChullType CHdeg(3); 
  for (it = V.begin(); it != V.end(); ++it)  CHdeg.insert(*it);

  V[0] = Point( RT(12), RT( 8), RT(-4), RT(2) );
  V[1] = Point( RT(12), RT( 8), RT(-4), RT(4) ); 
  V[2] = Point( RT(12), RT( 8), RT(-4), RT(8) );
  V[3] = Point( RT(12), RT( 8), RT(-4), RT(16) );
  V[4] = Point( RT(12), RT( 8), RT(-4), RT(1) );
  V[5] = Point( RT(12), RT( 8), RT(-4), RT(32) );
  V[6] = Point( RT(12), RT( 8), RT(-4), RT(2) );
  V[7] = Point( RT(12), RT( 8), RT(-4), RT(16) );

  ChullType CHdegdeg(3); 
  for (it = V.begin(); it != V.end(); ++it)  CHdegdeg.insert(*it);

  V[0] = Point( RT(12), RT( 8), RT(-4), RT(2) );
  V[1] = Point( RT(12), RT( 8), RT(-4), RT(2) ); 
  V[2] = Point( RT(12), RT( 8), RT(-4), RT(2) );
  V[3] = Point( RT(12), RT( 8), RT(-4), RT(2) );
  V[4] = Point( RT(12), RT( 8), RT(-4), RT(2) );
  V[5] = Point( RT(12), RT( 8), RT(-4), RT(2) );
  V[6] = Point( RT(12), RT( 8), RT(-4), RT(2) );
  V[7] = Point( RT(12), RT( 8), RT(-4), RT(2) );

  ChullType CHdegdegdeg(3); 
  for (it = V.begin(); it != V.end(); ++it)  CHdegdegdeg.insert(*it);

  return true;
}
#endif // _TEST_CLS_CHULL_TRAITS_3_C
