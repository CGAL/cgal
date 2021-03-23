// Copyright (c) 2020 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_ONE_RING_NEIGHBOR_QUERY_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_ONE_RING_NEIGHBOR_QUERY_H

#include <CGAL/license/Shape_detection.h>

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>

namespace CGAL {
namespace Shape_detection {
namespace Polyline {

  template<
  typename GeomTraits,
  typename InputRange>
  class One_ring_neighbor_query {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;

    One_ring_neighbor_query(
      const InputRange& input_range) :
    m_input_range(input_range) {

      CGAL_precondition(input_range.size() > 0);
    }

    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      neighbors.clear();
      const std::size_t n = m_input_range.size();
      CGAL_precondition(query_index < n);
      const std::size_t im = (query_index + n - 1) % n;
      const std::size_t ip = (query_index + 1) % n;
      neighbors.push_back(im);
      neighbors.push_back(ip);
    }

  private:
    const Input_range& m_input_range;
  };

} // namespace Polyline
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_ONE_RING_NEIGHBOR_QUERY_H
