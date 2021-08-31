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
// Author(s)     : Dmitry Anisimov, Gennadii Sytov
//

#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_WRAPPER_2_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_WRAPPER_2_H

#include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<typename GeomTraits>
  struct Segment_wrapper_2 {

    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2  = typename Traits::Direction_2;

    std::size_t index = std::size_t(-1); // regularize segments and contours
    bool is_used = false;

    Point_2 barycenter, ref_coords; // regularize segments
    FT orientation, length, a, b, c;
    Direction_2 direction;

    Segment_2 segment; // regularize contours
    std::size_t group = std::size_t(-1);
    bool is_valid_direction = false;

    void set_index(const std::size_t index_) {
      index = index_; is_used = true;
    }

    void set_length(
      const Segment_2& segment) {
      length = internal::segment_length_2(segment);
    }

    void set_barycenter(
      const Segment_2& segment) {
      barycenter = CGAL::midpoint(segment.source(), segment.target());
    }

    void set_direction(
      const Segment_2& segment) {
      auto v = segment.to_vector();
      direction = internal::direction_2(v);
    }

    void set_orientation(
      const Segment_2& segment) {
      set_direction(segment);
      orientation = internal::orientation_2(direction);
    }

    void set_abc(
      const Segment_2& segment) {
      set_barycenter(segment);
      set_orientation(segment);
      internal::line_coefficients_2(
        barycenter, direction, a, b, c);
    }

    void set_qp(
      const std::size_t index_,
      const Segment_2& segment) {
      set_index(index_);
      set_length(segment);
      set_abc(segment);
    }

    void set_ref_coords(
      const Point_2& frame_origin) {
      ref_coords = internal::transform_coordinates_2(
        barycenter, frame_origin, orientation);
    }
  };

  template<typename GeomTraits>
  class Wrap_segment_map {

    using Self = Wrap_segment_map<GeomTraits>;
    using key_type   = Segment_wrapper_2<GeomTraits>;
    using value_type = typename GeomTraits::Segment_2;
    using reference  = const value_type&;
    using category   = boost::readable_property_map_tag;

    reference operator[](const key_type& key) const {
      return key.segment;
    }
    friend reference get(const Self& self, const key_type& key) {
      return self[key];
    }
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_WRAPPER_2_H
