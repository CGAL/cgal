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

#include <CGAL/Visibility_2/visibility_utils.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/use.h>
#include <algorithm>
#include <iterator>
#include <string>

namespace CGAL {

/*! Computes the position of a watchtower of minimum vertical height above a
  1.5D terrain, from which the terrain is k-crossing visible. */
template <class Arrangement_2_, class RegularizationCategory = CGAL::Tag_true> class Watchtower_k_crossing_visibility
{

public:
  typedef Arrangement_2_ Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2 Geometry_traits_2;
  typedef typename Geometry_traits_2::Kernel Kernel;
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef CGAL::Polygon_2<Kernel> Polygon;

  typedef RegularizationCategory Regularization_category;

  Watchtower_k_crossing_visibility() : m_k(0) {}

  explicit Watchtower_k_crossing_visibility(unsigned const int k)
      : m_k(k) {}

  static std::string name() { return std::string("W_visibility_2"); }

  unsigned int k() const { return m_k; } // Number of allowed crossings

  void set_k(unsigned const int k) { m_k = k; }

  /*! Places a watchtower of minimum vertical height above the 1.5D terrain
    [first, beyond), an x-monotone polygonal chain given by its vertices in
    increasing x order, such that every point of the terrain is k-crossing
    visible from it. */
  template <class ForwardIterator>
  Point_2 compute_continuous_watchtower(ForwardIterator first, ForwardIterator beyond) const {
    // 1. construct the enclosing polygon T'
    const Polygon t_prime = construct_enclosing_polygon(first, beyond);

    // TODO
    CGAL_USE(t_prime);

    return Point_2();
  }

  /*! Places a watchtower of minimum vertical height above the 1.5D terrain
    [first, beyond), an x-monotone polygonal chain given by its vertices in
    increasing x order, such that every vertex of the terrain is k-crossing
    visible from it. */
  template <class ForwardIterator>
  Point_2 compute_discrete_watchtower(ForwardIterator first, ForwardIterator beyond) const {
    // TODO
    CGAL_USE(first);
    CGAL_USE(beyond);

    return Point_2();
  }

private:
  template <class ForwardIterator>
  Polygon construct_enclosing_polygon(ForwardIterator first, ForwardIterator beyond) const {
    CGAL_precondition(std::distance(first, beyond) >= 2);

    Polygon t_prime(first, beyond);

    FT y_max = first->y();
    FT x_right = first->x();
    for(ForwardIterator it = first; it != beyond; ++it) {
      y_max = (std::max)(y_max, it->y());
      x_right = it->x();
    }

    // TODO: h_T only has to lie above every point that the algorithm may ever
    // consider as a watchtower; pin down the bound this offset must satisfy.
    const FT y_top = y_max + FT(1);

    t_prime.push_back(Point_2(x_right, y_top));
    t_prime.push_back(Point_2(first->x(), y_top));

    CGAL_postcondition(t_prime.is_simple());
    CGAL_postcondition(t_prime.is_counterclockwise_oriented());

    return t_prime;
  }

  unsigned int m_k;
};

} // namespace CGAL

#endif // CGAL_WATCHTOWER_K_CROSSING_VISIBILITY_2_H
