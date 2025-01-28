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
// Author(s)     : Dmitry Anisimov, Simon Giraudot
//

#ifndef CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H

#include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Closed_contour_2.h>
#include <CGAL/Shape_regularization/internal/Open_contour_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  struct CLOSED { };
  struct OPEN { };

  template<
  typename ContourTag,
  typename ContourDirections,
  typename GeomTraits>
  class Contour_regularization_2 {

  public:
    using Contour_tag = ContourTag;
    using Contour_directions = ContourDirections;
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Regularization = typename std::conditional<
      std::is_same<ContourTag, CLOSED>::value,
      internal::Closed_contour_2<Contour_directions, Traits>,
      internal::Open_contour_2<Contour_directions, Traits> >::type;

    template<
    typename InputRange,
    typename PointMap,
    typename NamedParameters>
    Contour_regularization_2(
      const ContourDirections& directions,
      const InputRange& input_range,
      const PointMap point_map,
      const NamedParameters& np,
      const GeomTraits&) {

      CGAL_precondition(input_range.size() >= 2);
      const FT max_offset_2 = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::maximum_offset), FT(1) / FT(2));
      m_regularization = std::make_unique<Regularization>(
        directions, max_offset_2);
      m_regularization->initialize(input_range, point_map);
    }

    template<typename OutputIterator>
    OutputIterator regularize(OutputIterator contour) {
      return m_regularization->regularize(contour);
    }

  private:
    std::unique_ptr<Regularization> m_regularization;
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
