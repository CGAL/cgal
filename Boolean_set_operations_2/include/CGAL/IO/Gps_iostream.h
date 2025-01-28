// Copyright (c) 2009  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_GPS_IOSTREAM_H
#define CGAL_GPS_IOSTREAM_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/disable_warnings.h>

#include <iostream>
#include <list>
#include <algorithm>

#include <CGAL/basic.h>
#include <CGAL/General_polygon_set_2.h>

namespace CGAL {

/*! Inserter operator for general polygons sets.
 * Inserts a general polygon set into an output stream.
 * \param os the output stream.
 * \param pgn_set the general polygon set.
 * \return the output stream.
 */
template <typename GeomTraits_, typename Dcel_>
std::ostream&
operator<<(std::ostream& os,
           const CGAL::General_polygon_set_2<GeomTraits_, Dcel_>& pgn_set) {
  using Geometry_traits_2 = GeomTraits_;
  using Dcel = Dcel_;
  using Gps = CGAL::General_polygon_set_2<Geometry_traits_2, Dcel>;
  using Pwh_2 = typename Gps::Polygon_with_holes_2;
  using Pgn_with_holes_container = std::list<Pwh_2>;

  Pgn_with_holes_container pwhs;
  pgn_set.polygons_with_holes(std::back_inserter(pwhs));
  std::cout << pgn_set.number_of_polygons_with_holes() << std::endl;
  std::copy(pwhs.begin(), pwhs.end(), std::ostream_iterator<Pwh_2>(os, "\n"));
  return os;
}

/*! Extractor operator for general polygons sets.
 * Extracts a general polygon set from an input stream.
 * \param is the input stream.
 * \param pgn_set the general polygon set.
 * \return the input stream.
 */
template <typename GeomTraits_, typename Dcel_>
std::istream&
operator>>(std::istream& is,
           CGAL::General_polygon_set_2<GeomTraits_, Dcel_>& pgn_set) {
  using Geometry_traits_2 = GeomTraits_;
  using Dcel = Dcel_;
  using Gps = CGAL::General_polygon_set_2<Geometry_traits_2, Dcel>;
  using Pwh_2 = typename Gps::Polygon_with_holes_2;

  int n;
  is >> n;
  for (int i = 0; i < n; ++i) {
    Pwh_2 pwh;
    is >> pwh;
    pgn_set.insert(pwh);
  }
  return is;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
