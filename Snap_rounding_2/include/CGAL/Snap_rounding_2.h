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

  #ifndef CGAL_NO_DEPRECATED_CODE
/*!
\deprecated This function is deprecated since \cgal 6.2, use newer API of \ref snap_rounding_2_fct "CGAL::snap_rounding_2()" or `CGAL::hot_pixel_snap_rounding_2()`

\ingroup Snap_rounding_common

\tparam Traits must be a model of `HotPixelSnapRoundingTraits_2`.
\tparam InputIterator must be an iterator with value type `Traits::Segment_2`.
\tparam OutputContainer must be a container with a method `push_back(const OutputContainer::value_type& c)`,
where `OutputContainer::value_type` must be a container with a method `push_back(const Traits::Point_2& p)`

\param begin,end denote the iterator range
of the input segments.

\param output_container is a reference to a
container of the output polylines. Since a polyline is composed of a
sequence of points, a polyline is a container itself.

\param pixel_size denotes the pixel size `w`. The plane will be
tiled with square pixels of width `w` such that the origin is the
center of a pixel. `pixel_size` must have a positive value.

\param do_isr determines whether to apply ISR or SR.

\param int_output If set to true, the output coordinates are expressed in the integer grid (pixel indices).
Otherwise, they are given in the input coordinate system.

\param number_of_kd_trees The seventh parameter is briefly described later on this page; for a
detailed description see \cgalCite{cgal:hp-isr-02}.
*/
template<class Traits, class InputIterator, class OutputContainer>
CGAL_DEPRECATED void snap_rounding_2(InputIterator begin,
                                     InputIterator end,
                                     OutputContainer & output_container,
                                     typename Traits::NT pixel_size,
                                     bool do_isr = true,
                                     bool int_output = true,
                                     unsigned int number_of_kd_trees = 1)
{
  Snap_Rounding_2::internal::hot_pixel_snap_rounding_2<Traits>(begin, end, output_container, pixel_size, do_isr, int_output, number_of_kd_trees);
}
#endif

#if DOXYGEN_RUNNING
/**
* \ingroup Snap_rounding_common
*
* \brief subdivides and rounds a range of segments so that they are pairwise disjoint in their interiors.
* \anchor snap_rounding_2_fct
*
* By default, each polyline of the output corresponds to an input segment. Consequently, duplicate segments may appear in the output, for instance when multiple input segments collapse.
* When the parameter `output_unique_segments` is set to `true`, the polylines are decomposed into individual segments (represented as polylines with two points), and duplicates are removed.
*
* By default, this function rounds on double precision coordinates using `CGAL::vertical_slab_snap_rounding_2()`.
* Other rounding schemes or methods can be used by providing a `geom_traits` that is model of `VerticalSlabSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`.
*
* @tparam SegmentRange model of the concept `ConstRange` whose iterator is model of `ForwardIterator` and whose value_type is `geom_traits::Segment_2`, where the type of `geom_traits` is detailed by `np::geom_traits`.
* @tparam OutputContainer model of the concept `BackInsertionSequence` whose value type is itself a model of the concepts `DefaultConstructible` and `BackInsertionSequence` whose value type is `geom_traits::Point_2`
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* \param segments the input segment range
* \param out the output container
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" including the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{The traits class must respect the concept of `VerticalSlabSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `CGAL::Double_grid_snap_rounding_traits_2<Kernel>` where Kernel is deduced from the segment type, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
*   \cgalParamNBegin{output_unique_segments}
*     \cgalParamDescription{If set to true, the output polylines are unique pairs of distinct points represented a segment. As a result, the total number of output polylines may differ from the number of input segments.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{true}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template <class SegmentRange , class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
OutputContainer snap_rounding_2(const SegmentRange &segments,
                                OutputContainer    &out,
                                const NamedParameters &np = parameters::default_values());

/**
* \ingroup Snap_rounding_common
*
* \brief subdivides and rounds a range of polygons so that their boundary segments are pairwise disjoint in their interiors.
*
* If the input polygons are disjoint, the output polygons remain non-overlapping, although they may share vertices or edges.
* Each output polygon is free of self-intersections but may present pinched sections.
*
* By default, this function rounds on double precision coordinates using `CGAL::vertical_slab_snap_rounding_2()`.
* Other rounding schemes or methods can be used by providing a `geom_traits` that is model of `VerticalSlabSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`.
*
* @tparam PolygonRange model of a ConstRange whose iterator is model of ForwardIterator and whose value_type is model of `CGAL::Polygon_2`.
* @tparam OutputContainer model of the concept `BackInsertionSequence` whose value type is model of `CGAL::Polygon_2`.
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
*     \cgalParamDefault{an instance of `CGAL::Double_grid_snap_rounding_traits_2<Kernel>` where %Kernel is deduced from the point type of the polygons, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @warning a convex input polygon might no longer be convex after rounding.
*/
template <class PolygonRange, class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
void snap_rounding_2(const PolygonRange &polygons,
                     OutputContainer    &out,
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

template <class InputRange , class OutputContainer, class NamedParameters = parameters::Default_named_parameters>
void snap_rounding_2(const InputRange &inputs,
                     OutputContainer  &out,
                     const NamedParameters &np = parameters::default_values())
{
  using parameters::is_default_parameter;

  if constexpr(is_default_parameter<NamedParameters, internal_np::geom_traits_t>::value)
  {
    vertical_slab_snap_rounding_2(inputs, out, np);
  }
  else
  {
    using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t, NamedParameters, int>::type;
    if constexpr(internal::is_instance_of_VerticalSlabSRT_v<Traits>)
    {
      vertical_slab_snap_rounding_2(inputs, out, np);
    }
    else
    {
      hot_pixel_snap_rounding_2(inputs, out, np);
    }
  }
}
#endif

} //namespace CGAL

#endif
