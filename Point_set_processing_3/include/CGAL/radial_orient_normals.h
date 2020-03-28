// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Laurent Saboret

#ifndef CGAL_RADIAL_ORIENT_NORMALS_H
#define CGAL_RADIAL_ORIENT_NORMALS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Origin.h>
#include <CGAL/IO/trace.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <deque>
#include <math.h>

namespace CGAL {

/// \cond SKIP_IN_MANUAL

/// Radial orientation of the `[first, beyond)` range of points.
/// Normals are oriented towards exterior of the point set.
/// This very fast method is intended to convex objects.
///
/// This method modifies the order of input points so as to pack all oriented points first,
/// and returns an iterator over the first point with an unoriented normal (see erase-remove idiom).
/// For this reason it should not be called on sorted containers.
///
/// \pre normals must be unit vectors
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
///        It can be omitted if the value type of `ForwardIterator` is convertible to `Point_3<Kernel>`.
/// @tparam NormalPMap is a model of `ReadWritePropertyMap` with value type `Vector_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointPMap`.
///
/// @return iterator over the first point with an unoriented normal.

// This variant requires all parameters.
template <typename ForwardIterator,
          typename PointPMap,
          typename NormalPMap,
          typename Kernel
>
ForwardIterator
radial_orient_normals(
    ForwardIterator first,  ///< iterator over the first input point.
    ForwardIterator beyond, ///< past-the-end iterator over the input points.
    PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
    NormalPMap normal_pmap, ///< property map: value_type of ForwardIterator -> Vector_3.
    const Kernel& kernel) ///< geometric traits.
{
    CGAL_TRACE("Calls radial_orient_normals()\n");

    // Input points types
    typedef typename std::iterator_traits<ForwardIterator>::value_type Enriched_point;
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename boost::property_traits<NormalPMap>::value_type Vector;
    typedef typename boost::property_traits<PointPMap>::reference Point_ref;
    typedef typename boost::property_traits<NormalPMap>::reference Vector_ref;
    typedef typename Kernel::FT FT;

    // Precondition: at least one element in the container.
    CGAL_point_set_processing_precondition(first != beyond);

    // Find points barycenter.
    // Note: We should use CGAL::centroid() from PCA component.
    //       Unfortunately, it is not compatible with property maps.
    Vector sum = CGAL::NULL_VECTOR;
    int nb_points = 0;
    for (ForwardIterator it = first; it != beyond; it++)
    {
      Point_ref point = get(point_pmap, *it);
      sum = sum + (point - CGAL::ORIGIN);
      nb_points++;
    }
    Point barycenter = CGAL::ORIGIN + sum / (FT)nb_points;

    // Iterates over input points and orients normals towards exterior of the point set.
    // Copy points with robust normal orientation to oriented_points[], the others to unoriented_points[].
    std::deque<Enriched_point> oriented_points, unoriented_points;
    for (ForwardIterator it = first; it != beyond; it++)
    {
      Point_ref point = get(point_pmap, *it);

      // Radial vector towards exterior of the point set
      Vector vec1 = point - barycenter;

      // Point's normal
      Vector_ref vec2 = get(normal_pmap, *it);

      //         ->               ->
      // Orients vec2 parallel to vec1
      double dot = vec1 * vec2;
      if (dot < 0)
        vec2 = -vec2;

      put(normal_pmap, *it, vec2);

      // Is orientation robust?
      bool oriented = (std::abs(dot) > std::cos(80.*CGAL_PI/180.)); // robust iff angle < 80 degrees
      if (oriented)
        oriented_points.push_back(*it);
      else
        unoriented_points.push_back(*it);
    }

    // Replaces [first, beyond) range by the content of oriented_points[], then unoriented_points[].
    ForwardIterator first_unoriented_point =
      std::copy(oriented_points.begin(), oriented_points.end(), first);
    std::copy(unoriented_points.begin(), unoriented_points.end(), first_unoriented_point);

    CGAL_TRACE("  => %u normals are unoriented\n", unoriented_points.size());
    CGAL_TRACE("End of radial_orient_normals()\n");

    return first_unoriented_point;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator,
          typename PointPMap,
          typename NormalPMap
>
ForwardIterator
radial_orient_normals(
    ForwardIterator first,  ///< iterator over the first input point.
    ForwardIterator beyond, ///< past-the-end iterator over the input points.
    PointPMap point_pmap, ///< property map: value_type of ForwardIterator -> Point_3.
    NormalPMap normal_pmap) ///< property map: value_type of ForwardIterator -> Vector_3.
{
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    return radial_orient_normals(
      first,beyond,
      point_pmap, normal_pmap,
      Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
template <typename ForwardIterator,
          typename NormalPMap
>
ForwardIterator
radial_orient_normals(
    ForwardIterator first,  ///< iterator over the first input point.
    ForwardIterator beyond, ///< past-the-end iterator over the input points.
    NormalPMap normal_pmap) ///< property map: value_type of ForwardIterator -> Vector_3.
{
    return radial_orient_normals(
      first,beyond,
      make_identity_property_map(
      typename std::iterator_traits<ForwardIterator>::value_type()),
      normal_pmap);
}
/// @endcond

/// @endcond


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_RADIAL_ORIENT_NORMALS_H
