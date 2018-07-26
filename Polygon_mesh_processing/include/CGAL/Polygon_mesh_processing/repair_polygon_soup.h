// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_POLYGON_SOUP
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_POLYGON_SOUP

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Kernel_traits.h>

#include <algorithm>
#include <iterator>
#include <map>
#include <vector>
#include <utility>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template <typename PointRange, typename PolygonRange>
struct Polygon_types
{
  typedef typename PointRange::value_type                                           Point_3;
  typedef typename PolygonRange::value_type                                         Polygon_3;

  typedef typename std::iterator_traits<typename Polygon_3::iterator>::value_type   V_ID;
  typedef typename std::vector<Polygon_3>::size_type                                P_ID;
};

template <typename Stream, typename Polygon>
void print_polygon(Stream& out, const Polygon& polygon)
{
  const std::size_t polygon_size = polygon.size();
  out << "(" << polygon_size << ")";
  for(std::size_t i=0; i<polygon_size; ++i)
    out << " " << polygon[i];
  out << std::endl;
}

template <typename Traits, typename PointRange, typename Polygon>
bool simplify_polygon(PointRange& points,
                      Polygon& polygon,
                      const Traits& traits = Traits())
{
  const std::size_t ini_polygon_size = polygon.size();

  // Start at the last since if two points are identical, the second one gets removed.
  // By starting at 'last', we ensure that 'to_remove' is ordered from closest to .begin()
  // to closest to .end()
  std::size_t last = ini_polygon_size - 1, i = last;
  bool stop = false;
  std::vector<std::size_t> to_remove;

  do
  {
    std::size_t next_i = (i == last) ? 0 : i+1;
    stop = (next_i == last);

    while(polygon[i] == polygon[next_i] || // combinatorial equality
          traits.equal_3_object()(points[polygon[i]], points[polygon[next_i]])) // geometric equality
    {
      to_remove.push_back(next_i);
#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_DEBUG
      std::cout << "removing point: polygon[" << next_i << "] = " << polygon[next_i] << std::endl;
#endif
      next_i = (next_i == last) ? 0 : next_i+1;

      // Every duplicate in front of 'last' (circularly-speaking) has already been cleared
      if(next_i == last)
      {
        stop = true;
        break;
      }
    }

    i = next_i;
  }
  while(!stop);

  while(!to_remove.empty())
  {
    polygon.erase(polygon.begin() + to_remove.back());
    to_remove.pop_back();
  }

  const std::size_t removed_points_n = ini_polygon_size - polygon.size();
  return (removed_points_n != 0);
}

// \ingroup PMP_repairing_grp
//
// For each polygon of the soup, removes consecutive identical (in a geometric sense) points.
//
// \tparam PointRange a model of the concept `RandomAccessContainer`
// \tparam PolygonRange a model of the concept `SequenceContainer`
//                      whose value_type is itself a model of the concept `SequenceContainer`
//                      whose value_type is `std::size_t`.
// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
//
// \param points points of the soup of polygons.
// \param polygons a vector of polygons. Each element in the vector describes a polygon
//        using the indices of the points in `points`.
// \param np optional \ref pmp_namedparameters "Named Parameters" described below
//
// \cgalNamedParamsBegin
//    \cgalParamBegin{geom_traits} a geometric traits class instance.
//       The traits class must provide the nested functor `Equal_3`
//       to compare lexicographically two points a function `Equal_3 equal_3_object()`.
//   \cgalParamEnd
// \cgalNamedParamsEnd
//
template <typename Traits, typename PointRange, typename PolygonRange>
std::size_t simplify_polygons_in_polygon_soup(PointRange& points,
                                              PolygonRange& polygons,
                                              const Traits& traits = Traits())
{
  typedef typename Polygon_types<PointRange, PolygonRange>::P_ID        P_ID;
  typedef typename Polygon_types<PointRange, PolygonRange>::Polygon_3   Polygon_3;

  std::size_t simplified_polygons_n = 0;

  for(P_ID polygon_index=0, end=polygons.size(); polygon_index!=end; ++polygon_index)
  {
    Polygon_3& polygon = polygons[polygon_index];
    if(polygon.size() <= 1)
      continue;

    if(simplify_polygon(points, polygon, traits))
      ++simplified_polygons_n;
  }

  return simplified_polygons_n;
}

// \ingroup PMP_repairing_grp
//
// Splits "pinched" poygons, that is polygons for which a point appears more than once,
// into multiple non-pinched polygons.
//
// \tparam PointRange a model of the concept `RandomAccessContainer`
// \tparam PolygonRange a model of the concept `SequenceContainer`
//                      whose value_type is itself a model of the concept `SequenceContainer`
//                      whose value_type is `std::size_t`.
// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
//
// \param points points of the soup of polygons.
// \param polygons a vector of polygons. Each element in the vector describes a polygon
//        using the indices of the points in `points`.
// \param np optional \ref pmp_namedparameters "Named Parameters" described below
//
// \cgalNamedParamsBegin
//    \cgalParamBegin{geom_traits} a geometric traits class instance.
//       The traits class must provide the nested functor `Less_xyz_3`
//       to compare lexicographically two points a function `Less_xyz_3 less_xyz_3_object()`.
//   \cgalParamEnd
// \cgalNamedParamsEnd
//
template <typename Traits, typename PointRange, typename PolygonRange>
std::size_t split_pinched_polygons_in_polygon_soup(PointRange& points,
                                                   PolygonRange& polygons,
                                                   const Traits& traits = Traits())
{
  typedef typename Polygon_types<PointRange, PolygonRange>::P_ID        P_ID;
  typedef typename Polygon_types<PointRange, PolygonRange>::Point_3     Point_3;
  typedef typename Polygon_types<PointRange, PolygonRange>::Polygon_3   Polygon_3;

  typedef typename Traits::Less_xyz_3                                   Less_xyz_3;

  std::size_t new_polygons_n = 0;

  // It is important that polygons.size() is re-evaluated at each loop iteration
  // because new polygons are intentionally added at the back of 'polygons' to also be examined
  for(P_ID polygon_index=0; polygon_index < polygons.size(); ++polygon_index)
  {
    Polygon_3& polygon = polygons[polygon_index];
    const std::size_t ini_polygon_size = polygon.size();

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_DEBUG
    std::cout << "Input polygon: ";
    print_polygon(std::cout, polygon);
#endif

    if(ini_polygon_size <= 3)
      continue;

    CGAL_assertion(!simplify_polygon(points, polygon, traits)); // 'polygon' must not have duplicates

    typedef std::map<Point_3, std::size_t, Less_xyz_3>                   Unique_point_container;
    Unique_point_container unique_points(traits.less_xyz_3_object());

    for(std::size_t i=0, polygon_size = polygon.size(); i<polygon_size; ++i)
    {
      const Point_3& p = points[polygon[i]];

      std::pair<typename Unique_point_container::iterator, bool> is_insert_succesful =
        unique_points.insert(std::make_pair(p, i));

      if(!is_insert_succesful.second)
      {
        // We have already met that point, split the polygon into two smaller polygons
        std::size_t prev_id = is_insert_succesful.first->second;

        CGAL_assertion(prev_id < i-1);

        Polygon_3 split_polygon_1(polygon.begin() + prev_id, polygon.begin() + i);
        Polygon_3 split_polygon_2; // might be pinched too, but it'll be checked later

        split_polygon_2.insert(split_polygon_2.end(), polygon.begin(), polygon.begin() + prev_id);
        split_polygon_2.insert(split_polygon_2.end(), polygon.begin() + i, polygon.end());

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_DEBUG
        std::cout << "New polygons:" << std::endl;
        std::cout << "P1:";
        print_polygon(std::cout, split_polygon_1);

        std::cout << "P2:";
        print_polygon(std::cout, split_polygon_2);
#endif

        std::swap(polygon, split_polygon_1);
        polygons.push_back(split_polygon_2);

        ++new_polygons_n;
        break;
      }
    }
  }

  return new_polygons_n;
}

// \ingroup PMP_repairing_grp
//
// Removes polygons with fewer than 2 points from the soup.
//
// \tparam PointRange a model of the concept `Container`
// \tparam PolygonRange a model of the concept `SequenceContainer`
//                      whose value_type is itself a model of the concept `Container`
//                      whose value_type is `std::size_t`.
// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
//
// \param points points of the soup of polygons.
// \param polygons a vector of polygons. Each element in the vector describes a polygon
//        using the indices of the points in `points`.
//
template <typename PointRange, typename PolygonRange>
std::size_t remove_invalid_polygons_in_polygon_soup(PointRange& /*points*/,
                                                    PolygonRange& polygons)
{
  std::vector<std::size_t> to_remove;
  const std::size_t ini_polygons_size = polygons.size();
  for(std::size_t polygon_index=0; polygon_index!=ini_polygons_size; ++polygon_index)
  {
    if(polygons[polygon_index].size() <= 2)
      to_remove.push_back(polygon_index);
  }

  while(!to_remove.empty())
  {
    polygons.erase(polygons.begin() + to_remove.back());
    to_remove.pop_back();
  }

  return ini_polygons_size - polygons.size();
}

} // end namespace internal

template <typename PointRange, typename PolygonRange>
std::size_t remove_degenerate_polygons_in_polygon_soup(PointRange& points,
                                                       PolygonRange& polygons)
{
  return remove_degenerate_polygons_in_polygon_soup(points, polygons, CGAL::parameters::all_default());
}

/// \ingroup PMP_repairing_grp
///
/// Removes the isolated points from a polygon soup.
/// A point is considered <i>isolated</i> if it does not appear in any polygon of the soup.
///
/// \tparam PointRange a model of the concept `SequenceContainer`
/// \tparam PolygonRange a model of the concept `RandomAccessContainer`
///                      whose value_type is itself a model of the concept `RandomAccessContainer`
///                      whose value_type is `std::size_t`.
///
/// \param points points of the soup of polygons.
/// \param polygons a vector of polygons. Each element in the vector describes a polygon
///        using the indices of the points in `points`.
///
/// \returns the number of removed isolated points
///
template <typename PointRange, typename PolygonRange>
std::size_t remove_isolated_points_in_polygon_soup(PointRange& points,
                                                   PolygonRange& polygons)
{
  typedef typename internal::Polygon_types<PointRange, PolygonRange>::P_ID        P_ID;
  typedef typename internal::Polygon_types<PointRange, PolygonRange>::Polygon_3   Polygon_3;

  if(points.empty())
    return 0;

  // Go through all the polygons to find points that are never used.
  const std::size_t ini_points_size = points.size();
  std::vector<bool> visited(ini_points_size, false);
  std::vector<std::size_t> id_remapping(ini_points_size);
  for(std::size_t i=0; i<ini_points_size; ++i)
    id_remapping[i] = i;

  for(P_ID polygon_index=0, end=polygons.size(); polygon_index!=end; ++polygon_index)
  {
    const Polygon_3& polygon = polygons[polygon_index];

    for(std::size_t i=0, polygon_size = polygon.size(); i<polygon_size; ++i)
      visited[polygon[i]] = true;
  }

  // Move all the unused points to the end
  std::size_t swap_position = ini_points_size - 1;
  for(std::size_t i=0; i<ini_points_size; ++i)
  {
    if(!visited[i])
    {
#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_DEBUG
      std::cout << "points[" << i << "] = " << points[i] << " is isolated" << std::endl;
#endif
      std::swap(points[swap_position], points[i]);
      id_remapping[swap_position] = i;
      --swap_position;
    }
  }

  // Actually remove the unused points
  ++swap_position; // 'swap_position' points at the first element to remove
  const std::size_t removed_points_n = ini_points_size - swap_position;
  points.erase(points.begin() + swap_position, points.end());

  // Renumber the polygons
  for(P_ID polygon_index=0, end=polygons.size(); polygon_index!=end; ++polygon_index)
  {
    Polygon_3& polygon = polygons[polygon_index];
    for(std::size_t i=0, polygon_size = polygon.size(); i<polygon_size; ++i)
    {
      polygon[i] = id_remapping[polygon[i]];
      CGAL_postcondition(polygon[i] < points.size());
    }
  }

  return removed_points_n;
}

/// \ingroup PMP_repairing_grp
///
/// Cleans a given polygon soup by removing invalid elements. More precisely, this function
/// carries out the following tasks, in the same order as they are listed:
/// - simplification of the polygons to remove identical consecutive points;
/// - splitting of "pinched" polygons, that is polygons where a position appears more than once,
///   in multiple non-pinched polygons.
/// - removal of invalid polygons, that is polygons with fewer than 2 points;
/// - removal of isolated points, that is points that do not appear in any polygon of the soup,
///   using `CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup()`.
///
/// \tparam PointRange a model of the concept `SequenceContainer`
/// \tparam PolygonRange a model of the concept `SequenceContainer`
///                      whose value_type is itself a model of the concept `SequenceContainer`
///                      whose value_type is `std::size_t`.
/// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// \param points points of the soup of polygons.
/// \param polygons a vector of polygons. Each element in the vector describes a polygon
///        using the indices of the points in `points`.
/// \param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{geom_traits} a geometric traits class instance.
///       The traits class must provide the nested functors :
///         - `Less_xyz_3` to compare lexicographically two points
///         - `Equal_3` to check whether 2 points are identical
///       and, for each functor `Foo`, a function `Foo foo_object()`.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename PointRange, typename PolygonRange, typename NamedParameters>
void repair_polygon_soup(PointRange& points,
                         PolygonRange& polygons,
                         const NamedParameters& np)
{
  using boost::get_param;
  using boost::choose_param;

  typedef typename boost::lookup_named_param_def <
    internal_np::geom_traits_t,
    NamedParameters,
    typename CGAL::Kernel_traits<typename internal::Polygon_types<PointRange, PolygonRange>::Point_3 >::type
  > ::type Traits;

  Traits traits = choose_param(get_param(np, internal_np::geom_traits), Traits());

  internal::simplify_polygons_in_polygon_soup(points, polygons, traits);
  internal::split_pinched_polygons_in_polygon_soup(points, polygons, traits);
  internal::remove_invalid_polygons_in_polygon_soup(points, polygons);

  remove_isolated_points_in_polygon_soup(points, polygons);
}

template <typename PointRange, typename PolygonRange>
void repair_polygon_soup(PointRange& points,
                         PolygonRange& polygons)
{
  return repair_polygon_soup(points, polygons, CGAL::parameters::all_default());
}

} // end namespace Polygon_mesh_processing

} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_POLYGON_SOUP
