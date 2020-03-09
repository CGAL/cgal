// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_CLUSTER_POINT_SET_H
#define CGAL_CLUSTER_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/squared_distance_3.h>
#include <CGAL/Point_set_processing_3/internal/Neighbor_query.h>
#include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>
#include <CGAL/for_each.h>

#include <queue>

namespace CGAL
{


/// \cond SKIP_IN_MANUAL
namespace Point_set_processing_3
{

namespace internal
{

// Trick to both compile version with Emptyset_iterator and with
// user-provided OutputIterator. Many output iterators (such as
// `std::back_insert_iterator`) cannot be default constructed, which
// makes the mechanism `choose_param(get_param(...),Default())` fails.
template <typename NamedParameters, typename OutputIterator>
OutputIterator get_neighborhood (const NamedParameters& np, OutputIterator*)
{
  return CGAL::parameters::get_parameter(np, internal_np::neighborhood);
}

template <typename NamedParameters>
CGAL::Emptyset_iterator get_neighborhood (const NamedParameters&, CGAL::Emptyset_iterator*)
{
  return CGAL::Emptyset_iterator();
}

} // namespace internal

} // namespace Point_set_processing_3
/// \endcond

// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms
   Identifies simply connected clusters on a nearest neighbors graph.

   \tparam PointRange is a model of `Range`. The value type of its
   iterator is the key type of the named parameter `point_map`.
   \tparam ClusterMap is a model of `ReadWritePropertyMap` with value
   type `std::size_t`.

   \param points input point range.
   \param cluster_map maps each point to the index of the cluster it belongs to.
   \param k number of neighbors.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{callback} an instance of
      `std::function<bool(double)>`. It is called regularly when the
      algorithm is running: the current advancement (between 0. and
      1.) is passed as parameter. If it returns `true`, then the
      algorithm continues its execution normally; if it returns
      `false`, the algorithm is stopped and the number of already
      computed clusters is returned.\cgalParamEnd
     \cgalParamBegin{neighbor_radius} spherical neighborhood radius. If
     provided, the neighborhood of a query point is computed with a fixed spherical
     radius instead of a fixed number of neighbors. In that case, the parameter
     `k` is used as a limit on the number of points returned by each spherical
     query (to avoid overly large number of points in high density areas). If no
     limit is wanted, use `k=0`.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \return the number of clusters identified.
*/
template <typename PointRange, typename ClusterMap, typename NamedParameters>
std::size_t cluster_point_set (PointRange& points,
                               ClusterMap cluster_map,
                               unsigned int k,
                               const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  
  // basic geometric types
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
  typedef typename Point_set_processing_3::GetNeighborhood<PointRange, NamedParameters>::type Neighborhood;
  typedef typename GetSvdTraits<NamedParameters>::type SvdTraits;

  CGAL_static_assertion_msg(!(boost::is_same<SvdTraits,
                              typename GetSvdTraits<NamedParameters>::NoTraits>::value),
                            "Error: no SVD traits");

  PointMap point_map = choose_parameter(get_parameter(np, internal_np::point_map), PointMap());
  typename Kernel::FT neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius),
                                                         typename Kernel::FT(0));
  typename Kernel::FT factor = choose_parameter(get_parameter(np, internal_np::attraction_factor),
                                                typename Kernel::FT(2));

  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                               std::function<bool(double)>());

  double callback_factor = 1.;
  if (!std::is_same<Neighborhood,
      typename Point_set_processing_3::GetNeighborhood<PointRange, NamedParameters>::Empty>::value)
    callback_factor = 0.5;

  typedef typename Kernel::Point_3 Point;

  // types for K nearest neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  // Init cluster map with -1
  for (const value_type& p : points)
    put (cluster_map, p, -1);

  Neighbor_query neighbor_query (points, point_map);

  std::queue<iterator> todo;
  std::size_t nb_clusters = 0;

  // Flooding algorithm from each point
  std::size_t done = 0;
  std::size_t size = points.size();
  
  for (iterator it = points.begin(); it != points.end(); ++ it)
  {
    const value_type& p = *it;
    
    if (get (cluster_map, p) != -1)
      continue;

    todo.push (it);

    while (!todo.empty())
    {
      iterator current = todo.front();
      todo.pop();

      if (get (cluster_map, *current) != -1)
        continue;

      put (cluster_map, *current, nb_clusters);
      ++ done;
      
      if (callback && !callback (callback_factor * (done + 1) / double(size)))
        return (nb_clusters + 1);

      neighbor_query.get_iterators (get (point_map, *current), k, neighbor_radius,
                                    boost::make_function_output_iterator
                                    ([&](const iterator& it) { todo.push(it); }));

    }

    ++ nb_clusters;
  }

  if (!std::is_same<Neighborhood,
      typename Point_set_processing_3::GetNeighborhood<PointRange, NamedParameters>::Empty>::value)
  {
    Neighborhood neighborhood = Point_set_processing_3::internal::get_neighborhood(np, (Neighborhood*)(nullptr));
    k *= factor;
    neighbor_radius *= factor;

    std::vector<iterator> neighbors;
    std::vector<std::pair<std::size_t, std::size_t> > adjacencies;

    done = 0;
    for (const value_type& p : points)
    {
      std::size_t c0 = get (cluster_map, p);
      
      neighbors.clear();
      neighbor_query.get_iterators (get (point_map, p), k, neighbor_radius,
                                    std::back_inserter (neighbors));

      for (const iterator& it : neighbors)
      {
        std::size_t c1 = get (cluster_map, *it);
        if (c0 < c1)
          adjacencies.push_back (std::make_pair (c0, c1));
        else if (c0 > c1)
          adjacencies.push_back (std::make_pair (c1, c0));
        // else c0 == c1, ignore
      }

      ++ done;
      if (callback && !callback (callback_factor + callback_factor * (done + 1) / double(size)))
        return nb_clusters;
    }
    std::sort (adjacencies.begin(), adjacencies.end());
    auto last = std::unique (adjacencies.begin(), adjacencies.end());
    std::copy (adjacencies.begin(), last, neighborhood);
  }
  
  return nb_clusters;
}

/// \cond SKIP_IN_MANUAL
// overload with default NP
template <typename PointRange, typename ClusterMap>
std::size_t cluster_point_set (PointRange& points,
                               ClusterMap cluster_map,
                               unsigned int k)
{
  return cluster_point_set (points, cluster_map, k,
                            CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond

} // namespace CGAL


#endif // CGAL_CLUSTER_POINT_SET_H
