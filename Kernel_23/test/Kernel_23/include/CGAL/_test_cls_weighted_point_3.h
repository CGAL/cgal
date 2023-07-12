// Copyright (c) 2017
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL__TEST_CLS_WEIGHTED_POINT_3_H
#define CGAL__TEST_CLS_WEIGHTED_POINT_3_H

#include <CGAL/Origin.h>
#include <CGAL/Point_3.h>
#include <CGAL/use.h>
#include <CGAL/Weighted_point_3.h>

#include <cassert>
#include <iostream>

using CGAL::internal::use;

template <class R>
bool _test_cls_weighted_point_3(const R& )
{
  std::cout << "Testing class Weighted_point_3" ;

  typedef typename  R::RT    RT;
  typedef typename  R::FT    FT;

  RT n1(-35);
  RT n2( 50);
  RT n3( 20);
  RT n4(  5);
  FT iw = -1;

  CGAL::Point_3<R> p0(n1, n2, n3, n4);

  // constructions
  typename R::Weighted_point_3       iwp;
  CGAL::Weighted_point_3<R> wp0; // default
  CGAL::Weighted_point_3<R> wp1(CGAL::ORIGIN); // with origin
  CGAL::Weighted_point_3<R> wp2(p0); // with Point_3
  CGAL::Weighted_point_3<R> wp3(p0, iw); // with Point_3 and weight
  CGAL::Weighted_point_3<R> wp4(iwp); // with R::Weighted_point_3
  CGAL::Weighted_point_3<R> wp5(wp3); // with CGAL::Weighted_point_3< R >
  CGAL::Weighted_point_3<R> wp6(n1, n2, n3); // with coordinates
  use(wp0); use(wp4); use(wp5);

  // assignment
  wp1 = wp6;

  std::cout << ".";

  // accessors
  CGAL::Point_3<R> p1 = wp2.point();
  assert(p1 == p0);

  FT w = wp3.weight();
  assert(w == iw);

  // no need to test the other operations as they use Point_3 operations (which
  // are tested in _test_cls_point_3.h)

  std::cout << "done" << std::endl;
  return true;
}

#endif // CGAL__TEST_CLS_WEIGHTED_POINT_3_H
