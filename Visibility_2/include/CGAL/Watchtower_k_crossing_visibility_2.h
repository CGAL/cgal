// Copyright (c) 2026 Toronto Metropolitan Unversity (Canada).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s):  Yeganeh Bahoo Torudi <bahoo@torontomu.ca>
//             Teodor Cirstoiu <tcirstoiu@torontomu.ca>
//

#ifndef CGAL_WATCHTOWER_K_CROSSING_VISIBILITY_2_H
#define CGAL_WATCHTOWER_K_CROSSING_VISIBILITY_2_H

#include <CGAL/license/Visibility_2.h>
#include <CGAL/use.h>
#include <iterator>

namespace CGAL {

/*! Places a watchtower of minimum vertical height above the 1.5D terrain
  [first, beyond), an x-monotone polygonal chain given by its vertices in
  increasing x order, such that every point of the terrain is k-crossing
  visible from it. */
template <class InputIterator>
 std::iterator_traits<InputIterator>::value_type
continuous_watchtower_k_crossing_visibility_2(InputIterator first,
                                              InputIterator beyond,
                                              unsigned int k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Point_2;

  // TODO
  CGAL_USE(first);
  CGAL_USE(beyond);
  CGAL_USE(k);

  return Point_2();
}

/*! Places a watchtower of minimum vertical height above the 1.5D terrain
  [first, beyond), an x-monotone polygonal chain given by its vertices in
  increasing x order, such that every vertex of the terrain is k-crossing
  visible from it. */
template <class InputIterator>
  std::iterator_traits<InputIterator>::value_type
discrete_watchtower_k_crossing_visibility_2(InputIterator first,
                                            InputIterator beyond,
                                            unsigned int k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Point_2;

  // TODO
  CGAL_USE(first);
  CGAL_USE(beyond);
  CGAL_USE(k);

  return Point_2();
}

} // namespace CGAL

#endif // CGAL_WATCHTOWER_K_CROSSING_VISIBILITY_2_H
