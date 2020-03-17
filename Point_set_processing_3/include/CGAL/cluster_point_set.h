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

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

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
OutputIterator get_adjacencies (const NamedParameters& np, OutputIterator*)
{
  return CGAL::parameters::get_parameter(np, internal_np::adjacencies);
}

template <typename NamedParameters>
CGAL::Emptyset_iterator get_adjacencies (const NamedParameters&, CGAL::Emptyset_iterator*)
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
   Identifies connected components on a nearest neighbors graph built
   using a query sphere of fixed radius centered on each point.

   \tparam PointRange is a model of `Range`. The value type of its
   iterator is the key type of the named parameter `point_map`.
   \tparam ClusterMap is a model of `ReadWritePropertyMap` with value
   type `std::size_t`.

   \param points input point range.
   \param cluster_map maps each point to the index of the cluster it belongs to.
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
     \cgalParamBegin{neighbor_radius} spherical neighborhood
     radius. If no value is provided, the default value is 1% of the
     bounding box diagonal.\cgalParamEnd
     \cgalParamBegin{attraction_factor} used to compute adjacencies
     between clusters. Adjacencies are computed using a nearest
     neighbor graph built similarly to the one used for clustering,
     using `attraction_factor * neighbor_radius` as
     parameter. %Default value is `2`.\cgalParamEnd
     \cgalParamBegin{adjacencies} model of `OutputIterator` that
     accepts objects of type `std::pair<std::size_t,
     std::size_t>`. Each pair contains the indices of two adjacent
     clusters. If this parameter is not used, adjacencies are not
     computed at all.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \return the number of clusters identified.
*/
template <typename PointRange, typename ClusterMap, typename NamedParameters>
std::size_t cluster_point_set (PointRange& points,
                               ClusterMap cluster_map,
                               const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  
  // basic geometric types
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
  typedef typename Point_set_processing_3::GetAdjacencies<PointRange, NamedParameters>::type Adjacencies;
  typedef typename GetSvdTraits<NamedParameters>::type SvdTraits;

  CGAL_static_assertion_msg(!(boost::is_same<SvdTraits,
                              typename GetSvdTraits<NamedParameters>::NoTraits>::value),
                            "Error: no SVD traits");

  PointMap point_map = choose_parameter(get_parameter(np, internal_np::point_map), PointMap());
  typename Kernel::FT neighbor_radius = choose_parameter(get_parameter(np, internal_np::neighbor_radius),
                                                         typename Kernel::FT(-1));
  typename Kernel::FT factor = choose_parameter(get_parameter(np, internal_np::attraction_factor),
                                                typename Kernel::FT(2));

  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                               std::function<bool(double)>());

  double callback_factor = 1.;
  if (!std::is_same<Adjacencies,
      typename Point_set_processing_3::GetAdjacencies<PointRange, NamedParameters>::Empty>::value)
    callback_factor = 0.5;

  typedef typename Kernel::Point_3 Point;

  // types for K nearest neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // If no radius is given, init with 1% of bbox diagonal
  if (neighbor_radius < 0)
  {
    CGAL::Bbox_3 bbox = CGAL::bbox_3 (CGAL::make_transform_iterator_from_property_map (points.begin(), point_map),
                                      CGAL::make_transform_iterator_from_property_map (points.end(), point_map));
    
    neighbor_radius = 0.01 * CGAL::approximate_sqrt
      ((bbox.xmax() - bbox.xmin()) * (bbox.xmax() - bbox.xmin())
       + (bbox.ymax() - bbox.ymin()) * (bbox.ymax() - bbox.ymin()) 
       + (bbox.zmax() - bbox.zmin()) * (bbox.zmax() - bbox.zmin()));
  }

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

      neighbor_query.get_iterators (get (point_map, *current), 0, neighbor_radius,
                                    boost::make_function_output_iterator
                                    ([&](const iterator& it) { todo.push(it); }), false);

    }

    ++ nb_clusters;
  }

  if (!std::is_same<Adjacencies,
      typename Point_set_processing_3::GetAdjacencies<PointRange, NamedParameters>::Empty>::value)
  {
    Adjacencies adjacencies = Point_set_processing_3::internal::get_adjacencies(np, (Adjacencies*)(nullptr));
    neighbor_radius *= factor;

    std::vector<iterator> neighbors;
    std::vector<std::pair<std::size_t, std::size_t> > adj;

    done = 0;
    for (const value_type& p : points)
    {
      std::size_t c0 = get (cluster_map, p);
      
      neighbors.clear();
      neighbor_query.get_iterators (get (point_map, p), 0, neighbor_radius,
                                    std::back_inserter (neighbors), false);

      for (const iterator& it : neighbors)
      {
        std::size_t c1 = get (cluster_map, *it);
        if (c0 < c1)
          adj.push_back (std::make_pair (c0, c1));
        else if (c0 > c1)
          adj.push_back (std::make_pair (c1, c0));
        // else c0 == c1, ignore
      }

      ++ done;
      if (callback && !callback (callback_factor + callback_factor * (done + 1) / double(size)))
        return nb_clusters;
    }
    std::sort (adj.begin(), adj.end());
    auto last = std::unique (adj.begin(), adj.end());
    std::copy (adj.begin(), last, adjacencies);
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
