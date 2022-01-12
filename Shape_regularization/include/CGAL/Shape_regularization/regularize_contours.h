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

#ifndef CGAL_SHAPE_REGULARIZATION_REGULARIZE_CONTOURS_H
#define CGAL_SHAPE_REGULARIZATION_REGULARIZE_CONTOURS_H

/// \cond SKIP_IN_MANUAL
#include <CGAL/license/Shape_regularization.h>
/// \endcond

/**
* \ingroup PkgShapeRegularizationRef
* \file CGAL/Shape_regularization/regularize_contours.h
* This header includes all classes for regularizing contours.
* It also includes the corresponding free functions.
*/

// Boost includes.
/// \cond SKIP_IN_MANUAL
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
/// \endcond

// Internal includes.
/// \cond SKIP_IN_MANUAL
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Contour_regularization_2.h>
/// \endcond

#include <CGAL/Shape_regularization/Contours/Longest_direction_2.h>
#include <CGAL/Shape_regularization/Contours/Multiple_directions_2.h>
#include <CGAL/Shape_regularization/Contours/User_defined_directions_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Contours {

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief regularizes closed contours.

    Given a set of ordered 2D points connected by segments, which form a closed contour,
    this function enables to reinforce three types of regularities among consecutive edges of this contour:
    - *Parallelism*: contour edges, which are detected as near parallel, are made exactly parallel.
    - *Orthogonality*: contour edges, which are detected as near orthogonal, are made exactly orthogonal.
    - *Collinearity*: parallel contour edges, which are detected as near collinear, are made exactly collinear.

    The principal directions of the contour are provided via the concept `ContourDirections`.

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam ContDirections
    a model of `ContourDirections`

    \tparam OutIterator
    a model of `OutputIterator` that accepts points of type `GeomTraits::Point_2`

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \param input_range
    a const range of ordered points, which form a contour

    \param directions
    estimated contour directions towards which the contour edges are oriented

    \param contour
    an output iterator with points of the regularized contour

    \param np
    an optional sequence of \ref bgl_namedparameters "Named Parameters"
    among the ones listed below; this parameter can be omitted,
    the default values are then used

    \cgalNamedParamsBegin
      \cgalParamNBegin{maximum_offset}
        \cgalParamDescription{maximum allowed orthogonal distance between two parallel
          and consecutive contour edges such that they are considered to be collinear}
        \cgalParamType{`GeomTraits::FT`}
        \cgalParamDefault{0.5 unit length}
      \cgalParamNEnd
      \cgalParamNBegin{point_map}
        \cgalParamDescription{a property map that maps an item from `input_range`
        to `GeomTraits::Point_2`}
        \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the input
        range and value type is `GeomTraits::Point_2`}
        \cgalParamDefault{`CGAL::Identity_property_map`}
      \cgalParamNEnd
      \cgalParamNBegin{geom_traits}
        \cgalParamDescription{an instance of geometric traits class}
        \cgalParamType{a model of `Kernel`}
        \cgalParamDefault{a CGAL `Kernel` deduced from the point type,
        using `CGAL::Kernel_traits`}
      \cgalParamNEnd
    \cgalNamedParamsEnd

    \return an output iterator to the element in the destination range,
    one past the last contour vertex stored

    \pre input_range.size() >= 3
  */
  template<
  typename InputRange,
  typename ContDirections,
  typename OutIterator,
  typename NamedParameters = parameters::Default_named_parameters>
  OutIterator regularize_closed_contour(
    const InputRange& input_range,
    const ContDirections& directions,
    OutIterator contour,
    const NamedParameters& np = parameters::default_values()) {

    using PointMap = typename CGAL::GetPointMap<InputRange, NamedParameters>::type;
    using Point_2 = typename PointMap::value_type;
    using GeomTraits = typename CGAL::Kernel_traits<Point_2>::Kernel;

    const PointMap point_map = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::point_map), PointMap());
    const GeomTraits traits = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::geom_traits), GeomTraits());

    CGAL_precondition(input_range.size() >= 3);
    using Contour_regularizer =
    internal::Contour_regularization_2<
      internal::CLOSED, ContDirections, GeomTraits>;

    Contour_regularizer regularizer(
      directions, input_range, point_map, np, traits);
    return regularizer.regularize(contour);
  }

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief regularizes closed contours.

    This function regularizes a closed contour with respect to the longest
    edge of this contour.
  */
  template<
  typename InputRange,
  typename OutIterator>
  OutIterator regularize_closed_contour(
    const InputRange& input_range,
    OutIterator contour) {

    CGAL_precondition(input_range.size() >= 3);
    using Iterator_type = typename InputRange::const_iterator;
    using Point_2 = typename std::iterator_traits<Iterator_type>::value_type;
    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    using Contour_directions = Longest_direction_2<GeomTraits, InputRange>;
    GeomTraits traits;

    Contour_directions directions(input_range, true);
    return regularize_closed_contour(
      input_range, directions, contour, CGAL::parameters::geom_traits(traits));
  }

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief regularizes open contours.

    Given a set of ordered 2D points connected by segments, which form an open contour,
    this function enables to reinforce three types of regularities among consecutive edges of this contour:
    - *Parallelism*: contour edges, which are detected as near parallel, are made exactly parallel.
    - *Orthogonality*: contour edges, which are detected as near orthogonal, are made exactly orthogonal.
    - *Collinearity*: parallel contour edges, which are detected as near collinear, are made exactly collinear.

    The principal directions of the contour are provided via the concept `ContourDirections`.

    \tparam InputRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam ContDirections
    a model of `ContourDirections`

    \tparam OutIterator
    a model of `OutputIterator` that accepts points of type `GeomTraits::Point_2`

    \tparam NamedParameters
    a sequence of \ref bgl_namedparameters "Named Parameters"

    \param input_range
    a const range of ordered points, which form a contour

    \param directions
    estimated contour directions towards which the contour edges are oriented

    \param contour
    an output iterator with points of the regularized contour

    \param np
    an optional sequence of \ref bgl_namedparameters "Named Parameters"
    among the ones listed below; this parameter can be omitted,
    the default values are then used

    \cgalNamedParamsBegin
      \cgalParamNBegin{maximum_offset}
        \cgalParamDescription{maximum allowed orthogonal distance between two parallel
          and consecutive contour edges such that they are considered to be collinear}
        \cgalParamType{`GeomTraits::FT`}
        \cgalParamDefault{0.5 unit length}
      \cgalParamNEnd
      \cgalParamNBegin{point_map}
        \cgalParamDescription{a property map that maps an item from `input_range`
        to `GeomTraits::Point_2`}
        \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type of the input
        range and value type is `GeomTraits::Point_2`}
        \cgalParamDefault{`CGAL::Identity_property_map`}
      \cgalParamNEnd
      \cgalParamNBegin{geom_traits}
        \cgalParamDescription{an instance of geometric traits class}
        \cgalParamType{a model of `Kernel`}
        \cgalParamDefault{a CGAL `Kernel` deduced from the point type,
        using `CGAL::Kernel_traits`}
      \cgalParamNEnd
    \cgalNamedParamsEnd

    \return an output iterator to the element in the destination range,
    one past the last contour vertex stored

    \pre input_range.size() >= 2
  */
  template<
  typename InputRange,
  typename ContDirections,
  typename OutIterator,
  typename NamedParameters = parameters::Default_named_parameters>
  OutIterator regularize_open_contour(
    const InputRange& input_range,
    const ContDirections& directions,
    OutIterator contour,
    const NamedParameters& np = parameters::default_values()) {

    using PointMap = typename CGAL::GetPointMap<InputRange, NamedParameters>::type;
    using Point_2 = typename PointMap::value_type;
    using GeomTraits = typename CGAL::Kernel_traits<Point_2>::Kernel;

    const PointMap point_map = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::point_map), PointMap());
    const GeomTraits traits = parameters::choose_parameter(
      parameters::get_parameter(np, internal_np::geom_traits), GeomTraits());

    CGAL_precondition(input_range.size() >= 2);
    using Contour_regularizer =
    internal::Contour_regularization_2<
      internal::OPEN, ContDirections, GeomTraits>;

    Contour_regularizer regularizer(
      directions, input_range, point_map, np, traits);
    return regularizer.regularize(contour);
  }

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief regularizes open contours.

    This function regularizes an open contour with respect to the longest
    edge of this contour.
  */
  template<
  typename InputRange,
  typename OutIterator>
  OutIterator regularize_open_contour(
    const InputRange& input_range,
    OutIterator contour) {

    CGAL_precondition(input_range.size() >= 3);
    using Iterator_type = typename InputRange::const_iterator;
    using Point_2 = typename std::iterator_traits<Iterator_type>::value_type;
    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    using Contour_directions = Longest_direction_2<GeomTraits, InputRange>;
    GeomTraits traits;

    Contour_directions directions(input_range, false);
    return regularize_open_contour(
      input_range, directions, contour, CGAL::parameters::geom_traits(traits));
  }

} // namespace Contours
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_REGULARIZE_CONTOURS_H
