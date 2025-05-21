// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     : Susan Hert


#ifndef CGAL__TEST_FCT_SEGMENT_3_H
#define CGAL__TEST_FCT_SEGMENT_3_H

template <class R>
bool
_test_fct_segment_3(const R& )
{
 std::cout << "Testing functions Segment_3" ;

 typedef typename  R::RT          RT;

 typedef typename  R::Point_3     Point_3;
 typedef typename  R::Segment_3   Segment_3;

 Point_3 p1 ( RT(0), RT(0), RT(0), RT(1) );
 Point_3 p2 ( RT(1), RT(1), RT(1), RT(1) );


  Segment_3 l1(p1, p2);
  assert( CGAL::squared_distance(l1.source(), l1.target()) == CGAL::squared_length(l1) );

 std::cout << "done" << std::endl;
 return true;
}

#endif // CGAL__TEST_FCT_SEGMENT_3_H
