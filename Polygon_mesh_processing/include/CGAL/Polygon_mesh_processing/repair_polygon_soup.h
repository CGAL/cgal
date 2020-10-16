// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_POLYGON_SOUP
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_POLYGON_SOUP

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/iterator.h>
#include <CGAL/Kernel_traits.h>

#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>
#include <boost/range.hpp>

#include <algorithm>
#include <iterator>
#include <ios>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <unordered_map>
#include <unordered_set>

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
  #ifndef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
    #define CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
  #endif
#endif

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template <typename PointRange, typename PolygonRange>
struct Polygon_types
{
  typedef typename boost::range_value<PointRange>::type                             Point_3;
  typedef typename boost::range_value<PolygonRange>::type                           Polygon_3;

  typedef typename boost::range_iterator<Polygon_3>::type                           V_ID_iterator;
  typedef typename std::iterator_traits<V_ID_iterator>::value_type                  V_ID;
  typedef typename std::vector<Polygon_3>::size_type                                P_ID;
};

template <typename PointRange, typename PolygonRange, typename NamedParameters>
struct GetPolygonGeomTraits
{
  typedef typename internal_np::Lookup_named_param_def <
                     internal_np::geom_traits_t,
                     NamedParameters,
                     typename CGAL::Kernel_traits<
                       typename internal::Polygon_types<
                         PointRange, PolygonRange>::Point_3 >::type
                   > ::type                                                         type;
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

template <typename Traits, typename PointRange>
struct Vertex_ID_comparer
{
  Vertex_ID_comparer(const PointRange& points, const Traits& traits = Traits())
    : points(points), traits(traits)
  { }

  template <typename VID>
  bool operator()(const VID id_1, const VID id_2) const {
    return traits.less_xyz_3_object()(points[id_1], points[id_2]);
  }

private:
  const PointRange& points;
  const Traits& traits;
};

template <typename Traits, typename PointRange, typename Polygon>
bool polygon_has_unique_vertices(const PointRange& points,
                                 const Polygon& polygon,
                                 const Traits& traits = Traits())
{
  typedef Vertex_ID_comparer<Traits, PointRange>                                  Comparer;

  Comparer comp(points, traits);
  std::set<std::size_t, Comparer> unique_vertices(comp);
  unique_vertices.insert(polygon.begin(), polygon.end());

  return (unique_vertices.size() == polygon.size());
}

template <typename Traits, typename PointRange, typename Polygon>
bool simplify_polygon(PointRange& points,
                      Polygon& polygon,
                      const Traits& traits = Traits())
{
  const std::size_t ini_polygon_size = polygon.size();

  for(std::size_t i=0; i<polygon.size(); ++i)
  {
    const std::size_t s = polygon.size();
    if(s == 1)
      break;

    const std::size_t ni = (i + 1) % s;
    if(polygon[i] == polygon[ni] ||
       traits.equal_3_object()(points[polygon[i]], points[polygon[ni]]))
    {
#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
      std::cout << "Duplicate point: polygon[" << ni << "] = " << polygon[ni] << std::endl;
#endif
      polygon.erase(polygon.begin() + i--);
    }
  }

  const std::size_t removed_points_n = ini_polygon_size - polygon.size();
  return (removed_points_n != 0);
}

// \ingroup PMP_repairing_grp
//
// For each polygon of the soup, removes consecutive identical (in a geometric sense) points.
//
// \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
// \tparam PolygonRange a model of the concept `SequenceContainer`
//                      whose value_type is itself a model of the concept `SequenceContainer`
//                      whose value_type is `std::size_t`.
// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
//
// \param points points of the soup of polygons.
// \param polygons a vector of polygons. Each element in the vector describes a polygon
//        using the indices of the points in `points`.
// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
//
// \cgalNamedParamsBegin
//   \cgalParamNBegin{geom_traits}
//     \cgalParamDescription{an instance of a geometric traits class}
//     \cgalParamType{The traits class must provide the nested functor `Equal_3`
//                    to compare lexicographically two points a function `Equal_3 equal_3_object()`.}
//     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
//     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
//   \cgalParamNEnd
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

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
  if(simplified_polygons_n > 0)
    std::cout << "Cleaned consecutive duplicate vertices in " << simplified_polygons_n << " polygon(s)" << std::endl;
#endif

  return simplified_polygons_n;
}

// \ingroup PMP_repairing_grp
//
// splits "pinched" polygons, that is polygons for which a point appears more than once,
// into multiple non-pinched polygons.
//
// \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
// \tparam PolygonRange a model of the concept `SequenceContainer`
//                      whose value_type is itself a model of the concepts `SequenceContainer`
//                      and `Swappable` whose value_type is `std::size_t`.
// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
//
// \param points points of the soup of polygons.
// \param polygons a vector of polygons. Each element in the vector describes a polygon
//        using the indices of the points in `points`.
// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
//
// \cgalNamedParamsBegin
//   \cgalParamNBegin{geom_traits}
//     \cgalParamDescription{an instance of a geometric traits class}
//     \cgalParamType{The traits class must provide the nested functor `Less_xyz_3`
//                    to compare lexicographically two points a function `Less_xyz_3 less_xyz_3_object()`.}
//     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
//     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
//   \cgalParamNEnd
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

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
    std::cout << "Input polygon: ";
    internal::print_polygon(std::cout, polygon);
#endif

    if(ini_polygon_size <= 3)
      continue;

    // 'polygon' must not have consecutive duplicates
    CGAL_assertion(!simplify_polygon(points, polygon, traits));

    typedef std::map<Point_3, std::size_t, Less_xyz_3>                   Unique_point_container;
    Unique_point_container unique_points(traits.less_xyz_3_object());

    for(std::size_t i=0, polygon_size = polygon.size(); i<polygon_size; ++i)
    {
      const Point_3& p = points[polygon[i]];

      std::pair<typename Unique_point_container::iterator, bool> is_insert_successful =
        unique_points.insert(std::make_pair(p, i));

      if(!is_insert_successful.second)
      {
#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
        std::cout << "Pinched polygon: ";
        internal::print_polygon(std::cout, polygon);
#endif

        // We have already met that point, split the polygon into two smaller polygons
        std::size_t prev_id = is_insert_successful.first->second;
        CGAL_assertion(prev_id < i-1);

        Polygon_3 split_polygon_1(polygon.begin() + prev_id, polygon.begin() + i);
        CGAL_postcondition(internal::polygon_has_unique_vertices(points, split_polygon_1, traits));

        Polygon_3 split_polygon_2; // might be pinched too, but it'll be checked later
        split_polygon_2.insert(split_polygon_2.end(), polygon.begin(), polygon.begin() + prev_id);
        split_polygon_2.insert(split_polygon_2.end(), polygon.begin() + i, polygon.end());

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
        std::cout << "New polygons:" << std::endl;
        std::cout << "P1:";
        internal::print_polygon(std::cout, split_polygon_1);

        std::cout << "P2:";
        internal::print_polygon(std::cout, split_polygon_2);
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
// removes polygons with fewer than 2 points from the soup.
//
// \tparam PointRange a model of the concept `Container` whose value type is the point type.
// \tparam PolygonRange a model of the concept `SequenceContainer`
//                      whose value_type is itself a model of the concept `Container`
//                      whose value_type is `std::size_t`.
// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
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
    {
#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
      std::cout << "Invalid polygon:";
      print_polygon(std::cout, polygons[polygon_index]);
#endif
      to_remove.push_back(polygon_index);
    }
  }

  while(!to_remove.empty())
  {
    polygons.erase(polygons.begin() + to_remove.back());
    to_remove.pop_back();
  }

  const std::size_t removed_polygons_n = ini_polygons_size - polygons.size();

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
  if(removed_polygons_n > 0)
    std::cout << "Removed " << removed_polygons_n << " invalid polygon(s)" << std::endl;
#endif

  return removed_polygons_n;
}

} // end namespace internal

/// \ingroup PMP_repairing_grp
///
/// removes the isolated points from a polygon soup.
/// A point is considered <i>isolated</i> if it does not appear in any polygon of the soup.
///
/// \tparam PointRange a model of the concept `SequenceContainer` whose value type is the point type.
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
  std::size_t first_unused_pos = ini_points_size;
  for(std::size_t i=0; i<ini_points_size;)
  {
    if(i >= first_unused_pos)
      break;

    if(!visited[i])
    {
      std::size_t swap_position = first_unused_pos - 1;
      CGAL_assertion(swap_position < ini_points_size);

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
      std::cout << "points[" << i << "] = " << points[i] << " is isolated" << std::endl;
      std::cout << "  swapping it to pos: " << swap_position << std::endl;
#endif
      std::swap(points[swap_position], points[i]);

      // Swap manually because MSVC is unhappy with std::swap(v[i], v[j]) and a vector<bool>,
      // and Apple Clang is unhappy with v.swap(v[i], v[j])...
      const bool tmp = visited[swap_position];
      visited[swap_position] = visited[i];
      visited[i] = tmp;

      id_remapping[swap_position] = i;
      --first_unused_pos;
    }
    else
    {
      ++i;
    }
  }

  // Actually remove the unused points
  const std::size_t removed_points_n = ini_points_size - first_unused_pos;

  // Pointless to remap everything if nothing has changed, so early exit
  if(removed_points_n == 0)
    return removed_points_n;

  points.erase(points.begin() + first_unused_pos, points.end());

  // Renumber the polygons
  for(P_ID polygon_index=0, end=polygons.size(); polygon_index!=end; ++polygon_index)
  {
    Polygon_3& polygon = polygons[polygon_index];
    for(std::size_t i=0, polygon_size = polygon.size(); i<polygon_size; ++i)
    {
      polygon[i] = id_remapping[polygon[i]];
      CGAL_postcondition(static_cast<std::size_t>(polygon[i]) < points.size());
    }
  }

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
  if(removed_points_n > 0)
    std::cout << "Removed " << removed_points_n << " isolated point(s)" << std::endl;
#endif

  return removed_points_n;
}

/// \ingroup PMP_repairing_grp
///
/// merges the duplicate points in a polygon soup.
/// Note that the index of a point that is merged with another point will thus change
/// in all the polygons that the point appears in.
///
/// \tparam PointRange a model of the concepts `SequenceContainer` and `Swappable`
///                    whose value type is the point type.
/// \tparam PolygonRange a model of the concept `RandomAccessContainer`
///                      whose value_type is itself a model of the concept `RandomAccessContainer`
///                      whose value_type is `std::size_t`.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param points points of the soup of polygons.
/// \param polygons a vector of polygons. Each element in the vector describes a polygon
///        using the indices of the points in `points`.
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested functor `Less_xyz_3`
///                    to compare lexicographically two points a function `Less_xyz_3 less_xyz_3_object()`.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \returns the number of removed points
///
template <typename PointRange, typename PolygonRange, typename NamedParameters>
std::size_t merge_duplicate_points_in_polygon_soup(PointRange& points,
                                                   PolygonRange& polygons,
                                                   const NamedParameters& np)
{
  typedef typename internal::Polygon_types<PointRange, PolygonRange>::P_ID        P_ID;
  typedef typename internal::Polygon_types<PointRange, PolygonRange>::Point_3     Point_3;
  typedef typename internal::Polygon_types<PointRange, PolygonRange>::Polygon_3   Polygon_3;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename internal::GetPolygonGeomTraits<PointRange, PolygonRange, NamedParameters>::type Traits;
  Traits traits = choose_parameter<Traits>(get_parameter(np, internal_np::geom_traits));

  typedef typename Traits::Less_xyz_3                                             Less_xyz_3;

  const std::size_t ini_points_n = points.size();
  std::vector<std::size_t> point_index(ini_points_n, 0);

  typedef std::map<Point_3, std::size_t, Less_xyz_3>                              Unique_point_container;
  Unique_point_container point_to_id(traits.less_xyz_3_object());

  std::vector<Point_3> unique_points;
  unique_points.reserve(ini_points_n);

  for(std::size_t i=0; i<ini_points_n; ++i)
  {
    std::pair<typename Unique_point_container::iterator, bool> is_insert_successful =
      point_to_id.insert(std::make_pair(points[i], unique_points.size()));

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
    if(!is_insert_successful.second)
      std::cout << "points[" <<i << "] = " << points[i] << " was already encountered" << std::endl;
#endif
    std::size_t id = is_insert_successful.first->second;

    if(id == unique_points.size())
      unique_points.push_back(points[i]);
    point_index[i] = id;
  }

  if(unique_points.size() != ini_points_n)
  {
    for(P_ID polygon_index=0, end=polygons.size(); polygon_index!=end; ++polygon_index)
    {
      Polygon_3& polygon = polygons[polygon_index];
#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
      std::cout << "Input polygon: ";
      internal::print_polygon(std::cout, polygon);
#endif

      for(std::size_t i=0, polygon_size = polygon.size(); i<polygon_size; ++i)
        polygon[i] = point_index[polygon[i]];

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
      std::cout << "Output polygon: ";
      internal::print_polygon(std::cout, polygon);
#endif
    }

    std::swap(points, unique_points);
  }

  const std::size_t removed_points_n = ini_points_n - points.size();

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
  if(removed_points_n > 0)
    std::cout << "Removed (merged) " << removed_points_n << " duplicate points" << std::endl;
#endif

  return removed_points_n;
}

template <typename PointRange, typename PolygonRange>
std::size_t merge_duplicate_points_in_polygon_soup(PointRange& points,
                                                   PolygonRange& polygons)
{
  return merge_duplicate_points_in_polygon_soup(points, polygons, CGAL::parameters::all_default());
}

namespace internal {

// Find the position of the (arbitrarily chose) first point of the canonical point
// and whether we should order from left to right or the opposite
template <typename Traits, typename PointRange, typename Polygon>
void canonical_polygon_markers(const PointRange& points,
                               const Polygon& polygon,
                               std::size_t& first,
                               bool& reversed,
                               const Traits& traits = Traits())
{
  CGAL_precondition(polygon.size() > 0);
  CGAL_precondition(polygon_has_unique_vertices(points, polygon, traits));

  typedef typename boost::range_iterator<const Polygon>::type                     V_ID_iterator;
  typedef Vertex_ID_comparer<Traits, PointRange>                                  Vertex_comparer;

  // Find the bottom-left-front-most point that will be the first point of the polygon
  Vertex_comparer comp(points, traits);
  V_ID_iterator min_id = std::min_element(polygon.begin(), polygon.end(), comp);
  first = min_id - polygon.begin();

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
  std::cout << "first: " << first
            << " points[" << *min_id << "] = " << points[*min_id] << std::endl;
#endif

  // Decide arbitrarily whether we are reading from left to right or the opposite
  const std::size_t last = polygon.size() - 1;
  std::size_t pos_prev = (first == 0) ? last : first - 1;
  std::size_t pos_next = (first == last) ? 0 : first + 1;

  reversed = traits.less_xyz_3_object()(points[polygon[pos_prev]], points[polygon[pos_next]]);

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
  std::cout << "pos_prev: " << pos_prev
            << " points[" << polygon[pos_prev] << "] = " << points[polygon[pos_prev]] << std::endl;
  std::cout << "pos_next: " << pos_next
            << " points[" << polygon[pos_next] << "] = " << points[polygon[pos_next]] << std::endl;
  std::cout << "reversed: " << std::boolalpha << reversed << std::endl;
#endif
}

template <typename Polygon>
Polygon construct_canonical_polygon_with_markers(const Polygon& polygon,
                                                 const std::size_t first,
                                                 const bool reversed)
{
  const std::size_t polygon_size = polygon.size();
  Polygon canonical_polygon;

  if(reversed)
  {
    std::size_t rfirst = polygon_size - 1 - first;
    canonical_polygon.insert(canonical_polygon.end(), polygon.rbegin() + rfirst, polygon.rend());
    canonical_polygon.insert(canonical_polygon.end(), polygon.rbegin(), polygon.rbegin() + rfirst);
  }
  else
  {
    canonical_polygon.insert(canonical_polygon.end(), polygon.begin() + first, polygon.end());
    canonical_polygon.insert(canonical_polygon.end(), polygon.begin(), polygon.begin() + first);
  }

  CGAL_postcondition(canonical_polygon[0] == polygon[first]);
  CGAL_postcondition(canonical_polygon.size() == polygon_size);
  return canonical_polygon;
}

// 'reversed' indicates whether the canonical polygon has the same order as input polygon.
template <typename Traits, typename PointRange, typename Polygon>
Polygon construct_canonical_polygon(const PointRange& points,
                                    const Polygon& polygon,
                                    bool& reversed,
                                    const Traits& traits = Traits())
{
  if(polygon.size() < 2)
  {
    reversed = false;
    return polygon;
  }


#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
  std::cout << "Input polygon:";
  internal::print_polygon(std::cout, polygon);
#endif

  std::size_t first;
  canonical_polygon_markers(points, polygon, first, reversed, traits);
  Polygon canonical_polygon = construct_canonical_polygon_with_markers(polygon, first, reversed);

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
  std::cout << "Canonical polygon:";
  internal::print_polygon(std::cout, canonical_polygon);
#endif

  return canonical_polygon;
}

template <typename Traits, typename PointRange, typename Polygon>
Polygon construct_canonical_polygon(const PointRange& points,
                                    const Polygon& polygon,
                                    const Traits& traits = Traits())
{
  bool useless = false;
  return construct_canonical_polygon(points, polygon, useless, traits);
}

template <typename PointRange, typename PolygonRange, typename Traits>
struct Polygon_hash
{
  typedef std::size_t                                                             result_type;
  typedef typename internal::Polygon_types<PointRange, PolygonRange>::Polygon_3   Polygon_3;

  Polygon_hash(const PointRange& points, const PolygonRange& canonical_polygons, const Traits& traits)
    : points(points), canonical_polygons(canonical_polygons), traits(traits)
  { }

  template <typename Polygon_ID>
  result_type operator() (const Polygon_ID& polygon_index) const
  {
    const Polygon_3& canonical_polygon = canonical_polygons[polygon_index];

    std::size_t seed = 0;
    for(std::size_t i=0, end=canonical_polygon.size(); i<end; ++i)
      boost::hash_combine(seed, canonical_polygon[i]);

    return seed;
  }

private:
  const PointRange& points;
  const PolygonRange& canonical_polygons;
  const Traits& traits;
};

template <typename PointRange, typename PolygonRange, typename Reversed_markers, typename Traits>
struct Polygon_equality_tester
{
  typedef bool                                                                    result_type;
  typedef typename internal::Polygon_types<PointRange, PolygonRange>::Polygon_3   Polygon_3;

  Polygon_equality_tester(const PointRange& points,
                          const PolygonRange& canonical_polygons,
                          const Reversed_markers& reversed_markers,
                          const Traits& traits,
                          const bool same_orientation = false)
    : points(points),
      canonical_polygons(canonical_polygons),
      reversed_markers(reversed_markers),
      traits(traits),
      same_orientation(same_orientation)
  { }

  template <typename Polygon_ID>
  bool operator()(const Polygon_ID& polygon_index_1, const Polygon_ID& polygon_index_2) const
  {
    const Polygon_3& canonical_polygon_1 = canonical_polygons[polygon_index_1];
    const Polygon_3& canonical_polygon_2 = canonical_polygons[polygon_index_2];

    if(same_orientation &&
       reversed_markers[polygon_index_1] != reversed_markers[polygon_index_2])
      return false;

    return (canonical_polygon_1 == canonical_polygon_2);
  }

private:
  const PointRange& points;
  const PolygonRange& canonical_polygons;
  const Reversed_markers& reversed_markers;
  const Traits& traits;
  const bool same_orientation;
};

template <typename ValueType, typename OutputIterator>
struct Duplicate_collector
{
  void collect_duplicates(const ValueType& v1, const ValueType& v2)
  {
    std::vector<ValueType>& verts = collections[v1];
    if(verts.empty())
      verts.push_back(v1);
    verts.push_back(v2);
  }

  void dump(OutputIterator out)
  {
    typedef std::pair<const ValueType, std::vector<ValueType> > Pair_type;
    for(const Pair_type& p : collections)
      *out++ = p.second;
  }

  std::unordered_map<ValueType, std::vector<ValueType> > collections;
};

template <typename ValueType>
struct Duplicate_collector<ValueType, CGAL::Emptyset_iterator>
{
  void collect_duplicates(const ValueType&, const ValueType&) { }
  void dump(CGAL::Emptyset_iterator) { }
};

// \ingroup PMP_repairing_grp
//
// collects duplicate polygons in a polygon soup, that is polygons that share the same vertices in the same
// order.
//
// \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
// \tparam PolygonRange a model of the concept `RandomAccessContainer`
//                      whose value_type is itself a model of the concepts `RandomAccessContainer`
//                      and `ReversibleContainer` whose value_type is `std::size_t`.
// \tparam DuplicateOutputIterator a model of `OutputIterator` with value type
//                                 `std::vector<std::vector<std::size_t> >`.
//
// \param points points of the soup of polygons.
// \param polygons a vector of polygons. Each element in the vector describes a polygon
//        using the indices of the points in `points`.
// \param out the output iterator in which duplicate polygons are put. Each entry is a vector of
//            polygon ids `i0`, `i1`, etc. such that `polygons[i0] = polygons[i1] = ...`
// \param same_orientation whether two polygons should have the same orientation to be duplicates.
//
template <typename PointRange, typename PolygonRange, typename DuplicateOutputIterator, typename Traits>
DuplicateOutputIterator collect_duplicate_polygons(const PointRange& points,
                                                   const PolygonRange& polygons,
                                                   DuplicateOutputIterator out,
                                                   const Traits& traits = Traits(),
                                                   const bool same_orientation = false)
{
  typedef typename internal::Polygon_types<PointRange, PolygonRange>::P_ID        P_ID;

  typedef internal::Polygon_hash<PointRange, PolygonRange, Traits>                Hasher;
  typedef boost::dynamic_bitset<>                                                 Reversed_markers;
  typedef internal::Polygon_equality_tester<PointRange, PolygonRange,
                                            Reversed_markers, Traits>             Equality;
  typedef std::unordered_set<P_ID, Hasher, Equality>                      Unique_polygons;

  const std::size_t polygons_n = polygons.size();

  // We want the hash function to return the same value if the polygons are the same,
  // regardless of circular permutations and different orientations.
  PolygonRange canonical_polygons(polygons_n);
  Reversed_markers is_reversed(polygons_n, 0);
  for(P_ID polygon_index=0, end=polygons.size(); polygon_index!=end; ++polygon_index)
  {
    bool reversed;
    canonical_polygons[polygon_index] =
      internal::construct_canonical_polygon(points, polygons[polygon_index], reversed, traits);

    if(reversed)
      is_reversed.set(polygon_index);
  }

  Hasher hash(points, canonical_polygons, traits);
  Equality equal(points, canonical_polygons, is_reversed, traits, same_orientation);

  Unique_polygons unique_polygons(polygons_n /*bucket size*/, hash, equal);
  Duplicate_collector<P_ID, DuplicateOutputIterator> duplicates;

  for(P_ID polygon_index=0, end=polygons.size(); polygon_index!=end; ++polygon_index)
  {
    std::pair<typename Unique_polygons::iterator, bool> is_insert_successful =
      unique_polygons.insert(polygon_index);

    if(!is_insert_successful.second)
    {
      const P_ID other_polygon_id = *(is_insert_successful.first);
#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
      std::cout << "polygon: " << polygon_index << " is a duplicate of polygon: " << other_polygon_id << std::endl;
#endif
      duplicates.collect_duplicates(other_polygon_id, polygon_index);
    }
  }

  duplicates.dump(out);
  return out;
}

} // end namespace internal

/// \ingroup PMP_repairing_grp
///
/// merges the duplicate polygons in a polygon soup. Two polygons are duplicate if they share the same
/// vertices in the same order. Note that the first vertex of the polygon does not matter, that is
/// the triangle `0,1,2` is a duplicate of the triangle `2,0,1`.
///
/// \tparam PointRange a model of the concept `RandomAccessContainer` whose value type is the point type.
/// \tparam PolygonRange a model of the concept `SequenceContainer`
///                      whose value_type is itself a model of the concepts `RandomAccessContainer`
///                      and `ReversibleContainer` whose value_type is `std::size_t`.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param points points of the soup of polygons.
/// \param polygons a vector of polygons. Each element in the vector describes a polygon
///        using the indices of the points in `points`.
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested functor `Less_xyz_3`
///                    to compare lexicographically two points a function `Less_xyz_3 less_xyz_3_object()`.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{erase_all_duplicates}
///     \cgalParamDescription{Parameter to indicate, when multiple polygons are duplicates,
///                           whether all the duplicate polygons should be removed
///                           or if one (arbitrarily chosen) face should be kept.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`false`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{require_same_orientation}
///     \cgalParamDescription{Parameter to indicate if polygon orientation should be taken
///                           into account when determining whether two polygons are duplicates,
///                           that is, whether e.g. the triangles `0,1,2` and `0,2,1` are duplicates.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`false`}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \returns the number of removed polygons
///
template <typename PointRange, typename PolygonRange, typename NamedParameters>
std::size_t merge_duplicate_polygons_in_polygon_soup(const PointRange& points,
                                                     PolygonRange& polygons,
                                                     const NamedParameters& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename internal::Polygon_types<PointRange, PolygonRange>::P_ID                         P_ID;

  const bool erase_all_duplicates = choose_parameter(get_parameter(np, internal_np::erase_all_duplicates), false);
  const bool same_orientation = choose_parameter(get_parameter(np, internal_np::require_same_orientation), false);

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
  std::cout << "Only polygons with the same orientation are duplicates: " << std::boolalpha << same_orientation << std::endl;
  std::cout << "Erase all duplicate polygons: " << std::boolalpha << erase_all_duplicates << std::endl;
#endif

  typedef typename internal::GetPolygonGeomTraits<PointRange, PolygonRange, NamedParameters>::type Traits;
  Traits traits = choose_parameter<Traits>(get_parameter(np, internal_np::geom_traits));

  std::vector<std::vector<P_ID> > all_duplicate_polygons;
  internal::collect_duplicate_polygons(points, polygons, std::back_inserter(all_duplicate_polygons), traits, same_orientation);

  if(all_duplicate_polygons.empty())
    return 0;

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
  std::cout << all_duplicate_polygons.size() << " duplicate(s)" << std::endl;
#endif

  // Move all polygons that will be removed to the end of container
  const std::size_t init_polygons_n = polygons.size();
  std::size_t swap_position = init_polygons_n - 1;

  std::vector<bool> treated(init_polygons_n, false);

  // PID_to_pos is to go from a polygon ID to its position in the polygons vector, and pos_to_PID
  // is to move the other way
  std::vector<std::size_t> PID_to_pos(init_polygons_n);
  std::vector<std::size_t> pos_to_PID(init_polygons_n);
  for(std::size_t i=0, ps=polygons.size(); i<ps; ++i)
  {
    PID_to_pos[i] = i;
    pos_to_PID[i] = i;
  }

  while(!all_duplicate_polygons.empty())
  {
    const std::vector<P_ID>& duplicate_polygons = all_duplicate_polygons.back();
    CGAL_assertion(duplicate_polygons.size() >= 2);

    std::size_t i = erase_all_duplicates ? 0 : 1;
    for(; i<duplicate_polygons.size(); ++i)
    {
      const P_ID polygon_to_remove_id = duplicate_polygons[i];
      if(treated[polygon_to_remove_id])
        continue;

      const P_ID polygon_to_remove_pos = PID_to_pos[polygon_to_remove_id];
      CGAL_assertion(swap_position < init_polygons_n);
      const P_ID polygon_at_swap_position_id = pos_to_PID[swap_position];

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE_PP
      std::cout << "Removing duplicate, PID: " << polygon_to_remove_id << " at position: " << polygon_to_remove_pos << std::endl;
      std::cout << "  swap position: " << swap_position << ", position of PID: " << polygon_at_swap_position_id << std::endl;
#endif

      // Need to keep track of who goes where
      PID_to_pos[polygon_at_swap_position_id] = polygon_to_remove_pos;
      PID_to_pos[polygon_to_remove_id] = swap_position;
      pos_to_PID[polygon_to_remove_pos] = polygon_at_swap_position_id;
      pos_to_PID[swap_position] = polygon_to_remove_id;

      CGAL_assertion(polygon_to_remove_pos <= swap_position);
      std::swap(polygons[swap_position], polygons[polygon_to_remove_pos]);
      --swap_position;

      treated[polygon_to_remove_id] = true;
    }

    all_duplicate_polygons.pop_back();
  }

  ++swap_position; // so that it points to the first removed polygon
  const std::size_t removed_polygons_n = init_polygons_n - swap_position;

  typename PolygonRange::iterator first = polygons.begin();
  std::advance(first, swap_position);
  polygons.erase(first, polygons.end());

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
  std::cout << "Removed " << removed_polygons_n << " duplicate polygon(s)" << std::endl;
  std::cout << polygons.size() << " polygon(s) left" << std::endl;
#endif

  return removed_polygons_n;
}

template <typename PointRange, typename PolygonRange>
std::size_t merge_duplicate_polygons_in_polygon_soup(PointRange& points,
                                                     PolygonRange& polygons)
{
  return merge_duplicate_polygons_in_polygon_soup(points, polygons, CGAL::parameters::all_default());
}

/// \ingroup PMP_repairing_grp
///
/// cleans a given polygon soup through various repairing operations. More precisely, this function
/// carries out the following tasks, in the same order as they are listed:
/// - merging of duplicate points, using the function
///   `CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup()`;
/// - simplification of polygons to remove geometrically identical consecutive vertices;
/// - splitting of "pinched" polygons, that is polygons in which a geometric position appears more than once.
///   The splitting process results in multiple non-pinched polygons;
/// - removal of invalid polygons, that is polygons with fewer than 2 vertices;
/// - removal of duplicate polygons, using the function
///   `CGAL::Polygon_mesh_processing::merge_duplicate_polygons_in_polygon_soup()`;
/// - removal of isolated points,
///   using the function `CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup()`.
///
/// Note that the point and polygon containers will be modified by the repairing operations,
/// and thus the indexation of the polygons will also be changed.
///
/// \tparam PointRange a model of the concepts `SequenceContainer` and `Swappable`
///                    and whose value type is the point type.
/// \tparam PolygonRange a model of the concept `SequenceContainer`.
///                      whose value_type is itself a model of the concepts `SequenceContainer`,
///                      `Swappable`, and `ReversibleContainer` whose value_type is `std::size_t`.
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param points points of the soup of polygons.
/// \param polygons a vector of polygons. Each element in the vector describes a polygon
///        using the indices of the points in `points`.
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested functors `Less_xyz_3` and `Equal_3`
///                    to respectivelycompare lexicographically two points and to check if 2 points
///                    are identical. For each functor `Foo`, a function `Foo foo_object()` must be provided.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{erase_all_duplicates}
///     \cgalParamDescription{Parameter to indicate, when multiple polygons are duplicates,
///                           whether all the duplicate polygons should be removed
///                           or if one (arbitrarily chosen) face should be kept.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`false`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{require_same_orientation}
///     \cgalParamDescription{Parameter to indicate if polygon orientation should be taken
///                           into account when determining whether two polygons are duplicates,
///                           that is, whether e.g. the triangles `0,1,2` and `0,2,1` are duplicates.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`false`}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
template <typename PointRange, typename PolygonRange, typename NamedParameters>
void repair_polygon_soup(PointRange& points,
                         PolygonRange& polygons,
                         const NamedParameters& np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename internal::GetPolygonGeomTraits<PointRange, PolygonRange, NamedParameters>::type Traits;
  Traits traits = choose_parameter<Traits>(get_parameter(np, internal_np::geom_traits));

#ifdef CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
  std::cout << "Repairing soup with " << points.size() << " points and " << polygons.size() << " polygons" << std::endl;
#endif

  merge_duplicate_points_in_polygon_soup(points, polygons, np);
  internal::simplify_polygons_in_polygon_soup(points, polygons, traits);
  internal::split_pinched_polygons_in_polygon_soup(points, polygons, traits);
  internal::remove_invalid_polygons_in_polygon_soup(points, polygons);
  merge_duplicate_polygons_in_polygon_soup(points, polygons, np);
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
