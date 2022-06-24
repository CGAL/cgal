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

  /*!
    \ingroup PkgShapeDetectionRGOnPolyline

    \brief Direct neighbors connectivity in a polyline.

    This class returns indices of the previous and next vertex
    in a polyline given as `InputRange.`

    \tparam GeomTraits
    a model of `Kernel`

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \cgalModels `NeighborQuery`
  */
  template<
  typename GeomTraits,
  typename InputRange>
  class One_ring_neighbor_query {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    /// \endcond

    /// Item type.
    using Item = typename InputRange::const_iterator;
    using Region = std::vector<Item>;

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param input_range
      an instance of `InputRange` with polyline vertices

      \pre `input_range.size() > 0`
    */
    One_ring_neighbor_query(
      const InputRange& input_range) :
    m_input_range(input_range) {

      CGAL_precondition(input_range.size() > 0);
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief implements `NeighborQuery::operator()()`.

      This operator retrieves the previous and next vertex in a polyline
      with respect to the vertex `query`.
      The `Items` are returned in `neighbors`.

      \param query
      `Item` of the query point

      \param neighbors
      `Items` of vertices, which are direct neighbors of the query point

      \pre `query_index < input_range.size()`
    */
    void operator()(
      const Item query,
      std::vector<Item>& neighbors) const {

      neighbors.clear();

      Item before;
      if (query == m_input_range.begin())
        before = m_input_range.end() - 1;
      else
        before = query - 1;

      Item after = query + 1;
      if (after == m_input_range.end())
        after = m_input_range.begin();

      neighbors.push_back(before);
      neighbors.push_back(after);
    }

    /// @}

  private:
    const Input_range& m_input_range;
  };

} // namespace Polyline
} // namespace Shape_detection
} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_POLYLINE_ONE_RING_NEIGHBOR_QUERY_H
