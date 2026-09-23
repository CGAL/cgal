// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_MESH_PROCESSING_HILBERT_SORT_H
#define CGAL_POLYGON_MESH_PROCESSING_HILBERT_SORT_H
#include <CGAL/license/Polygon_mesh_processing/core.h>

#include <CGAL/hilbert_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>

#include <numeric>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

template <typename Polygon, typename PointRange>
struct Point_of_polygon_property_map
{
  typedef Point_of_polygon_property_map<Polygon,PointRange> Self;

  typedef Polygon key_type;
  typedef typename PointRange::value_type value_type;
  typedef const value_type& reference;
  typedef boost::lvalue_property_map_tag category;

  Point_of_polygon_property_map(const PointRange& points)
  : points(points)
  {}

  const value_type& operator[](const key_type& array) const { return points[array[0]]; }

  friend reference get(const Self& pm, const key_type& array) { return pm.points[array[0]]; }

  private:
  const PointRange& points;
};



template <class RandomAccessContainer, class Permutation>
void permute_in_place(RandomAccessContainer& a, Permutation& permutation) {
  const std::size_t n = a.size();

  for (std::size_t i = 0; i < n; ++i) {
    while (permutation[i] != i) {
      std::size_t j = permutation[i];
      std::swap(a[i], a[j]);
      std::swap(permutation[i], permutation[j]);
    }
  }
}

/**
 * \ingroup PMP_misc_grp
 *
 * sorts the points and polygons of a polygon soup using the Hilbert sort.
 * \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
 * \tparam PolygonRange a model of the concepts `RandomAccessContainer`
 *                      whose `value_type` is itself a model of the concepts `RandomAccessContainer`
 *                      awhose `value_type` is an unsigned integer type
 *                      convertible to `std::size_t`
 * \param points points of the soup of polygons
 * \param polygons each element in the range describes a polygon using the indices of the vertices.
 */
template <typename PointRange, typename PolygonRange>
void hilbert_sort_polygon_soup(PointRange& points,
                               PolygonRange& polygons)
{
  using Point_3 = typename PointRange::value_type;
  using Polygon = typename PolygonRange::value_type;
  using Kernel = typename Kernel_traits<Point_3>::type;
  using Point_traits_3 =  Spatial_sort_traits_adapter_3<Kernel,
                                                        typename Pointer_property_map<Point_3>::type >;
  using Polygon_traits_3 = Spatial_sort_traits_adapter_3<Kernel,
  Point_of_polygon_property_map<Polygon, PointRange>>;

  std::vector<std::size_t> indices(points.size());

  std::iota(indices.begin(), indices.end(), 0);
  hilbert_sort( indices.begin(),
                indices.end(),
                Point_traits_3(make_property_map(points)));

  const std::vector<std::size_t> permutation(indices);

  permute_in_place(points, indices);

  for (Polygon& polygon : polygons) {
    for(int i = 0; i < 3; ++i) {
      polygon[i] = permutation[polygon[i]];
    }
  }

  Point_of_polygon_property_map<Polygon, PointRange> polygon_map(points);
  hilbert_sort(polygons.begin(),
               polygons.end(),
               Polygon_traits_3(polygon_map));
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_HILBERT_SORT_H
