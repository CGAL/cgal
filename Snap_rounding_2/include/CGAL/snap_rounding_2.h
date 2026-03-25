// Copyright (c) 2025 Geometry Factory.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Léo Valque

#ifndef CGAL_SNAP_ROUNDING_2_H
#define CGAL_SNAP_ROUNDING_2_H

#include <CGAL/license/Snap_rounding_2.h>

#include<CGAL/vertical_slab_snap_rounding_2.h>
#include<CGAL/hot_pixel_snap_rounding_2.h>

namespace CGAL {

#if DOXYGEN_RUNNING
/**
* \ingroup PkgSnapRounding2Ref
*
* \brief subdivides and rounds a range of segments so that they are pairwise disjoint in their interiors.
*
* The output is a range of polylines, where each polyline corresponds to an input segment.
*
* calls the function `CGAL::vertical_slab_snap_rounding_2()` or `CGAL::hot_pixel_snap_rounding_2()` depending if
* `geom_traits` is model of `VerticalSlabSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`.
*
* @tparam SegmentRange a range of whose value type is model of `Kernel::Segment_2`
* @tparam OutputPolylineIterator model of OutputIterator holding `Polyline`. `Polyline` must be a type that provides a `push_back(Point_2)` function.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param segments the input segment range
* \param out the output inserter
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" including the one listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `VerticalSlabSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `CGAL::Double_grid_snap_rounding_traits_2`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template <class SegmentRange , class OutputPolylineIterator, class NamedParameters = parameters::Default_named_parameters>
OutputPolylineIterator snap_rounding_2(const SegmentRange &segments,
                                       OutputPolylineIterator   out,
                                       const NamedParameters &np = parameters::default_values());

/**
* \ingroup PkgSnapRounding2Ref
*
* \brief subdivides and rounds a range of segments so that they are pairwise disjoint in their interiors.
*
* The output is a range of segments.
*
* calls the function `CGAL::vertical_slab_snap_rounding_2()` or `CGAL::hot_pixel_snap_rounding_2()` depending if
* the parameter `geom_traits` is model of `VerticalSlabSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`.
*
* @tparam SegmentRange a range whose value type is model of `Kernel::Segment_2`
* @tparam OutputSegmentIterator model of OutputIterator holding `Kernel::Segment_2`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param segments the input segment range
* \param out the output inserter
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" including the one listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `VerticalSlabSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `CGAL::Double_grid_snap_rounding_traits_2`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template <class SegmentRange , class OutputSegmentIterator, class NamedParameters = parameters::Default_named_parameters>
OutputSegmentIterator snap_rounding_2(const SegmentRange& segments,
                                      OutputSegmentIterator    out,
                                      const NamedParameters &np = parameters::default_values());

/**
* \ingroup PkgSnapRounding2Ref
*
* \brief subdivides and rounds a range of polygons so that their boundary segments are pairwise disjoint in their interiors.
*
* If the input polygons are disjoint, the output polygons remain non-overlapping, although they may share vertices or edges.
* Each output polygon is free of self-intersections but may present pinched sections.
*
* calls the function `CGAL::vertical_slab_snap_rounding_2()` or `CGAL::hot_pixel_snap_rounding_2()` depending if
* `geom_traits` is model of `VerticalSlabSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`.
*
* @tparam PolygonRange a range of `CGAL::Polygon_2`
* @tparam OutputPolygonIterator model of OutputIterator holding `CGAL::Polygon_2`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param polygons the range of input polygons
* \param out the output inserter
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" including the one listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `VerticalSlabSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `CGAL::Double_grid_snap_rounding_traits_2`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
* @warning a convex input polygon might no longer be convex after rounding.
*/
template <class PolygonRange, class OutputPolygonIterator, class NamedParameters = parameters::Default_named_parameters>
OutputPolygonIterator snap_rounding_2(PolygonRange  &polygons,
                                      OutputPolygonIterator out,
                                      const NamedParameters &np = parameters::default_values());

#else

namespace internal{
// If Construct_rounded_point_2 is defined by the Traits, it is model of VerticalSlabSRT
template <typename, typename = void>
struct is_instance_of_VerticalSlabSRT : std::false_type {};

template <typename T>
struct is_instance_of_VerticalSlabSRT<T, std::void_t<typename T::Construct_rounded_point_2>> : std::true_type {};

template <typename T>
inline constexpr bool is_instance_of_VerticalSlabSRT_v = is_instance_of_VerticalSlabSRT<T>::value;
}

template <class InputRange , class OutputIterator, class NamedParameters = parameters::Default_named_parameters>
OutputIterator snap_rounding_2(const InputRange &inputs,
                               OutputIterator   out,
                               const NamedParameters &np = parameters::default_values())
{
  using parameters::is_default_parameter;

  if constexpr(is_default_parameter<NamedParameters, internal_np::geom_traits_t>::value)
  {
    return vertical_slab_snap_rounding_2(inputs, out, np);
  }
  else
  {
    using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t, NamedParameters, int>::type;
    if constexpr(internal::is_instance_of_VerticalSlabSRT_v<Traits>)
    {
      return vertical_slab_snap_rounding_2(inputs, out, np);
    }
    else
    {
      return hot_pixel_snap_rounding_2(inputs, out, np);
    }
  }
}
#endif

} //namespace CGAL

#endif
