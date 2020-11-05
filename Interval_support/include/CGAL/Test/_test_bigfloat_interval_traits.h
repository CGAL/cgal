// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany),
// INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>
//
//    \brief provides test functions for the \c Bigfloat_Interval concept.

#include <CGAL/basic.h>
#include <CGAL/Testsuite/use.h>

#include <cstddef>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>

#include <cassert>
#include <CGAL/tags.h>
#include <CGAL/use.h>

#include <CGAL/Bigfloat_interval_traits.h>

#ifndef CGAL_TEST_BIGFLOAT_INTERVAL_TRAITS_H
#define CGAL_TEST_BIGFLOAT_INTERVAL_TRAITS_H

namespace CGAL {

template <class Bigfloat_interval_>
void test_bigfloat_interval_traits() {

  typedef typename CGAL::Bigfloat_interval_traits<Bigfloat_interval_>::Self BFIT;
  typedef typename BFIT::Type                                               BFI;
  CGAL_USE_TYPE(typename BFIT::Bound);

  typedef typename BFIT::Is_bigfloat_interval Is_bigfloat_interval;
  CGAL_USE_TYPE(Is_bigfloat_interval);
  // using CGAL::Tag_true;
  CGAL_static_assertion(( ::boost::is_same< Is_bigfloat_interval, CGAL::Tag_true>::value));

  const typename BFIT::Construct construct = typename BFIT::Construct();
  const typename BFIT::Set_precision set_precision = typename BFIT::Set_precision();
  const typename BFIT::Get_precision get_precision = typename BFIT::Get_precision();
  const typename BFIT::Relative_precision relative_precision
    = typename BFIT::Relative_precision();

  CGAL_USE(construct);
  CGAL_USE(set_precision);
  CGAL_USE(get_precision);
  CGAL_USE(relative_precision);

  CGAL::set_precision(BFI(),15);
  assert(CGAL::get_precision(BFI())    == 15);
  assert(CGAL::set_precision(BFI(),23) == 15);
  assert(CGAL::set_precision(BFI(),70) == 23);

  //TODO: define what get_significant_bits should do and test is. Better name ?
  // just a compile check
  assert(CGAL::relative_precision(CGAL::hull(BFI(2),BFI(1)))==0);
  assert(CGAL::relative_precision(CGAL::hull(BFI(3),BFI(1)))==-1);
  assert(CGAL::relative_precision(CGAL::hull(BFI(4),BFI(1)))==-2);
  assert(CGAL::relative_precision(CGAL::hull(BFI(5),BFI(1)))==-2);

  assert(CGAL::relative_precision(CGAL::hull(BFI(2),BFI(1)))==0);
  assert(CGAL::relative_precision(CGAL::hull(BFI(3),BFI(2)))==1);
  assert(CGAL::relative_precision(CGAL::hull(BFI(5),BFI(4)))==2);
  assert(CGAL::relative_precision(CGAL::hull(BFI(9),BFI(8)))==3);
  assert(CGAL::relative_precision(CGAL::hull(BFI(8),BFI(7)))==2);
}

} //namespace CGAL

#endif // CGAL_TEST_REAL_COMPARABLE_H
