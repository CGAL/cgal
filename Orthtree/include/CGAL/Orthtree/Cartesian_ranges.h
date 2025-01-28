// Copyright (c) 2020  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_ORTHTREE_CARTESIAN_RANGE_H
#define CGAL_ORTHTREE_CARTESIAN_RANGE_H

#include <CGAL/license/Orthtree.h>

#include <CGAL/Iterator_range.h>
#include <boost/iterator/zip_iterator.hpp>

namespace CGAL
{

namespace Orthtrees
{

namespace internal
{

template <typename Traits>
struct Cartesian_ranges
{
  using Point = typename Traits::Point_d;
  using Cartesian_const_iterator = typename Traits::Cartesian_const_iterator_d;

  using Range_single = CGAL::Iterator_range<Cartesian_const_iterator>;

  Range_single operator() (const Point& p) const
  {
    return CGAL::make_range (p.cartesian_begin(), p.cartesian_end());
  }

  using Range_pair
  = CGAL::Iterator_range
    <boost::zip_iterator
     <boost::tuple<Cartesian_const_iterator, Cartesian_const_iterator> > >;

  Range_pair operator() (const Point& a, const Point& b) const
  {
    return CGAL::make_range
      (boost::make_zip_iterator
       (boost::make_tuple (a.cartesian_begin(), b.cartesian_begin())),
        boost::make_zip_iterator
       (boost::make_tuple (a.cartesian_end(), b.cartesian_end())));
  }

};

} // namespace internal

} // namespace Orthtrees

} // namespace CGAL


#endif // CGAL_ORTHTREE_CARTESIAN_RANGE_H
