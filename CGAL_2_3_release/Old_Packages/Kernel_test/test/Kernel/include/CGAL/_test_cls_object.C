// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : _test_cls_object.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL__TEST_CLS_OBJECT_C
#define CGAL__TEST_CLS_OBJECT_C

#include <CGAL/_test_cls_object.h>

template <class R>
bool
_test_cls_object(const R&)
{
  typedef typename  R::RT   RT;
  std::cout << "testing class object" ;

  CGAL::Object o1;
  assert( o1.is_empty() );
  CGAL::Object o2;
  assert( o2.is_empty() );
  CGAL::Object o3;
  assert( o3.is_empty() );
  CGAL::Point_2<R>   p21( RT(1), RT(1) );
  CGAL::Point_2<R>   p22;
  CGAL::Point_2<R>   p23;
  CGAL::Point_2<R>   p24( RT(4), RT(4) );
  CGAL::Point_3<R>   p31;
  CGAL::Point_3<R>   p32;
  CGAL::Line_2<R>    l21( p21, p24);
  CGAL::Line_2<R>    l22;
  CGAL::Line_3<R>    l31;
  CGAL::Line_3<R>    l32;

  std::cout << '.';

  o1 = o2;
  o2 = CGAL::make_object(p21);
  assert( o1.is_empty() );
  o1 = CGAL::make_object(p31);
  assert(   CGAL::assign( p22, o2 ) );
  assert( p22 == p21 );
  assert( ! CGAL::assign( p32, o2 ) );
  assert( ! CGAL::assign( l22, o2 ) );
  o3 = o2;
  assert(   CGAL::assign( p23, o3 ) );
  assert( p23 == p21 );
  assert( ! CGAL::assign( p32, o2 ) );
  assert( ! CGAL::assign( l32, o2 ) );
  assert(   CGAL::assign( p32, o1 ) );
  assert( ! CGAL::assign( p21, o1 ) );
  assert( p21 == p23 );

  std::cout << '.';

  o2 = CGAL::make_object(l21);
  assert(   CGAL::assign( l22, o2 ) );
  assert( l22 == l21 );
  assert( ! CGAL::assign( l32, o2 ) );
  assert( ! CGAL::assign( p22, o2 ) );
  o3 = CGAL::make_object(l31);
  assert(   CGAL::assign( l32, o3 ) );
  assert( ! CGAL::assign( p32, o3 ) );

  std::cout << "done" << std::endl;
  return true;
}
#endif // CGAL__TEST_CLS_OBJECT_C
