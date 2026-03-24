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

#include<CGAL/Vertical_slabs_snap_rounding_2.h>
#include<CGAL/Hot_pixel_snap_rounding_2.h>

namespace CGAL {

#ifndef CGAL_NO_DEPRECATED_CODE
/*!
\deprecated use newer API of `CGAL::snap_rounding_2()` or `CGAL::hot_pixel_snap_rounding_2()`

\ingroup PkgSnapRounding2Ref

\tparam Traits must be a model of `SnapRoundingTraits_2`.
\tparam InputIterator must be an iterator with value type `Traits::Segment_2`.
\tparam OutputContainer must be a container with a method `push_back(const OutputContainer::value_type& c)`,
where `OutputContainer::value_type` must be a container with a method `push_back(const Traits::Point_2& p)`

\param begin,end The first two parameters denote the iterator range
of the input segments.

\param output_container is a reference to a
container of the output polylines. Since a polyline is composed of a
sequence of points, a polyline is a container itself.

\param do_isr The fifth parameter determines whether to apply ISR or SR.

\param pixel_size The fourth parameter denotes the pixel size `w`. The plane will be
tiled with square pixels of width `w` such that the origin is the
center of a pixel. `pixel_size` must have a positive value.

\param int_output The sixth parameter denotes the output
representation. If the value of the sixth parameter is `true`
then the centers of pixels constitute the integer grid, and hence the
vertices of the output polylines will be integers. For example, the
coordinates of the center of the pixel to the right of the pixel
containing the origin will be `(1,0)` regardless of the pixel width.
If the value of the sixth parameter is `false`, then the centers
of hot pixels (and hence the vertices of the output polylines) will
bear their original coordinates, which may not necessarily be
integers. In the latter case, the coordinates of the center of the
pixel to the right of the pixel containing the origin, for example,
will be `(w,0)`.

\param number_of_kd_trees The seventh parameter is briefly described later on this page; for a
detailed description see \cgalCite{cgal:hp-isr-02}.

Snap Rounding (SR, for short) is a well-known method for converting
arbitrary-precision arrangements of segments into a fixed-precision
representation \cgalCite{gght-srlse-97}, \cgalCite{gm-rad-98}, \cgalCite{h-psifp-99}. In
the study of robust geometric computing, it can be classified as a
finite precision approximation technique. Iterated Snap Rounding (ISR,
for short) is a modification of SR in which each vertex is at least
half-the-width-of-a-pixel away from any non-incident edge
\cgalCite{cgal:hp-isr-02}. This package supports both methods. Algorithmic
details and experimental results are given in \cgalCite{cgal:hp-isr-02}.

Given a finite collection \f$ \Sc\f$ of segments in the plane, the
arrangement of \f$ \Sc\f$ denoted \f$ \Ac(\Sc)\f$ is the subdivision of the plane
into vertices, edges, and faces induced by \f$ \Sc\f$. A <I>vertex</I> of the arrangement is either a segment endpoint or
the intersection of two segments. Given an arrangement of segments
whose vertices are represented with arbitrary-precision coordinates,
SR proceeds as follows. We tile the plane
with a grid of unit squares, <I>pixels</I>, each centered at a point
with integer coordinates. A pixel is <I>hot</I> if it contains a
vertex of the arrangement. Each vertex of the arrangement is replaced
by the center of the hot pixel containing it and each edge \f$ e\f$ is
replaced by the polygonal chain through the centers of the hot pixels
met by \f$ e\f$, in the same order as they are met by \f$ e\f$.

In a snap-rounded arrangement, the distance between a vertex and
a non-incident edge can be extremely small compared with the width of a
pixel in the grid used for rounding. ISR is a modification of SR which
makes a vertex and a non-incident edge well separated (the distance
between each is at least half-the-width-of-a-pixel). However, the
guaranteed quality of the approximation in ISR degrades. For more
details on ISR see \cgalCite{cgal:hp-isr-02}.

The traits used here must support (arbitrary-precision) rational number type as
this is a basic requirement of SR.



\cgalHeading{About the Number of kd-Trees}

A basic query used in the algorithm is to report the hot pixels of
size \f$ w\f$ that a certain segment \f$ s\f$ intersects. An alternative way to
do the same is to query the hot pixels' centers contained in a
Minkowski sum of \f$ s\f$ with a pixel of width \f$ w\f$ centered at the origin;
we denote this Minkowski sum by \f$ M(s)\f$. Since efficiently implementing
this kind of query is difficult, we use an orthogonal range-search
structure instead. We query with the bounding box \f$ B(M(s))\f$ of \f$ M(s)\f$
in a two-dimensional kd-tree which stores the centers of hot
pixels. Since \f$ B(M(s))\f$ in general is larger than \f$ M(s)\f$, we still
need to filter out the hot pixels which do not intersect \f$ s\f$.

While this approach is easy to implement with \cgal, it may incur
considerable overhead since the area of \f$ B(M(s))\f$ may be much larger
than the area of \f$ M(s)\f$, possibly resulting in many redundant hot pixels
to filter out. Our heuristic solution, which we describe next, is to
use a cluster of kd-trees rather than just one. The cluster includes
several kd-trees, each has the plane, and hence the centers of hot
pixels, rotated by a different angle in the first quadrant of the
plane; for our purpose, a rotation by angles outside this quadrant
is symmetric to a rotation by an angle in the first quadrant.

Given a parameter \f$ c\f$, the angles of rotation are \f$ (i - 1)
\frac{\pi}{2c}, i=1,\ldots,c\f$, and we construct a kd-tree
corresponding to each of these angles. Then for a query segment \f$ s\f$,
we choose the kd-tree for which the area of \f$ B(M(s))\f$ is the smallest,
in order to (potentially) get less hot pixels to filter out. Since
constructing many kd-trees may be costly, our algorithm avoids
building a kd-tree which it expects to be queried a relatively small
number of times (we estimate this number in advance). How many
kd-trees should be used? It is difficult to provide a simple
answer for that. There are inputs for which the time to build more
than one kd-tree is far greater than the time saved by having to
filter out less hot pixels (sparse arrangements demonstrate this
behavior), and there are inputs which benefit from using several
kd-trees. Thus, the user can control the number of kd-trees
with the parameter  `number_of_kd_trees`. Typically, but not
always, one kd-tree (the default) is sufficient.
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
  internal::hot_pixel_snap_rounding_2<Traits>(begin, end, output_container, pixel_size, do_isr, int_output, number_of_kd_trees);
}
#endif

#if DOXYGEN_RUNNING
/**
* \ingroup PkgSnapRounding2Ref
*
* \brief subdivides and rounds a range of segments so that they are pairwise disjoint in their interiors.
*
* The output is a range of polylines, where each polyline corresponds to an input segment.
*
* calls the function `CGAL::vertical_slabs_snap_rounding_2()` or `CGAL::hot_pixel_snap_rounding_2()` depending if
* `geom_traits` is model of `VerticalSlabsSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`.
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
*     \cgalParamType{The traits class must respect the concept of `VerticalSlabsSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `Double_snap_rounding_traits_2`}
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
* calls the function `CGAL::vertical_slabs_snap_rounding_2()` or `CGAL::hot_pixel_snap_rounding_2()` depending if
* the parameter `geom_traits` is model of `VerticalSlabsSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`.
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
*     \cgalParamType{The traits class must respect the concept of `VerticalSlabsSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`}
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
* Each polygon in the output are free of self intersections but may present pinched sections or/and common vertices or segments with
* other polygons.
*
* calls the function `CGAL::vertical_slabs_snap_rounding_2()` or `CGAL::hot_pixel_snap_rounding_2()` depending if
* `geom_traits` is model of `VerticalSlabsSnapRoundingTraits_2` or `HotPixelSnapRoundingTraits_2`.
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
*     \cgalParamType{The traits class must respect the concept of `VerticalSlabsSnapRoundingTraits_2`}
*     \cgalParamDefault{an instance of `Double_snap_rounding_traits_2`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
* @warning an input convex polygon might no longer be convex after rounding.
*/
template <class PolygonRange, class OutputPolygonIterator, class NamedParameters = parameters::Default_named_parameters>
OutputPolygonIterator snap_rounding_2(PolygonRange  &polygons,
                                      OutputPolygonIterator out,
                                      const NamedParameters &np = parameters::default_values());

#else

namespace internal{
// If Construct_rounded_point_2 is defined by the Traits, it is model of VerticalSlabsSRT
template <typename, typename = void>
struct is_instance_of_VerticalSlabsSRT : std::false_type {};

template <typename T>
struct is_instance_of_VerticalSlabsSRT<T, std::void_t<typename T::Construct_rounded_point_2>> : std::true_type {};

template <typename T>
inline constexpr bool is_instance_of_VerticalSlabsSRT_v = is_instance_of_VerticalSlabsSRT<T>::value;
}

template <class InputRange , class OutputIterator, class NamedParameters = parameters::Default_named_parameters>
OutputIterator snap_rounding_2(const InputRange &inputs,
                               OutputIterator   out,
                               const NamedParameters &np = parameters::default_values())
{
  using parameters::is_default_parameter;

  if constexpr(is_default_parameter<NamedParameters, internal_np::geom_traits_t>::value)
  {
    return vertical_slabs_snap_rounding_2(inputs, out, np);
  }
  else
  {
    using Traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t, NamedParameters, int>::type;
    if constexpr(internal::is_instance_of_VerticalSlabsSRT_v<Traits>)
    {
      return vertical_slabs_snap_rounding_2(inputs, out, np);
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
