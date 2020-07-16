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
#include <CGAL/Point_set_processing_3/internal/bbox_diagonal.h>
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
   Identifies connected components on a nearest neighbor graph built
   using a query sphere of fixed radius centered on each point.

   \tparam PointRange is a model of `Range`. The value type of its
   iterator is the key type of the named parameter `point_map`.
   \tparam ClusterMap is a model of `ReadWritePropertyMap` with value
   type `std::size_t`.

   \param points input point range.
   \param cluster_map maps each point to the index of the cluster it belongs to.
   \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`}
     \cgalParamNEnd

     \cgalParamNBegin{callback}
       \cgalParamDescription{a mechanism to get feedback on the advancement of the algorithm
                             while it's running and to interrupt it if needed}
       \cgalParamType{an instance of `std::function<bool(double)>`.}
       \cgalParamDefault{unused}
       \cgalParamExtra{It is called regularly when the
                       algorithm is running: the current advancement (between 0. and
                       1.) is passed as parameter. If it returns `true`, then the
                       algorithm continues its execution normally; if it returns
                       `false`, the algorithm is stopped and the number of already
                       computed clusters is returned.}
       \cgalParamExtra{The callback will be copied and therefore needs to be lightweight.}
     \cgalParamNEnd

     \cgalParamNBegin{neighbor_radius}
       \cgalParamDescription{the spherical neighborhood radius}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`1` percent of the bounding box diagonal}
     \cgalParamNEnd

     \cgalParamNBegin{attraction_factor}
       \cgalParamDescription{used to compute adjacencies between clusters.
                             Adjacencies are computed using a nearest neighbor graph built similarly
                             to the one used for clustering, using `attraction_factor * neighbor_radius` as
                             parameter.}
       \cgalParamType{floating scalar value}
       \cgalParamDefault{`2`}
     \cgalParamNEnd

     \cgalParamNBegin{adjacencies}
       \cgalParamDescription{an output iterator used to output pairs containing the indices of two adjacent clusters.}
       \cgalParamType{a model of `OutputIterator` that accepts objects of type `std::pair<std::size_t, std::size_t>`}
       \cgalParamDefault{`CGAL::Emptyset_iterator`}
       \cgalParamExtra{If this parameter is not used, adjacencies are not computed at all.}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
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
  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
  typedef typename Point_set_processing_3::GetAdjacencies<PointRange, NamedParameters>::type Adjacencies;

  CGAL_static_assertion_msg(!(boost::is_same<typename GetSvdTraits<NamedParameters>::type,
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

  // types for K nearest neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // If no radius is given, init with 1% of bbox diagonal
  if (neighbor_radius < 0)
    neighbor_radius = 0.01 * Point_set_processing_3::internal::bbox_diagonal (points, point_map);

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

    if (int(get (cluster_map, p)) != -1)
      continue;

    todo.push (it);

    while (!todo.empty())
    {
      iterator current = todo.front();
      todo.pop();

      if (int(get (cluster_map, *current)) != -1)
        continue;

      put (cluster_map, *current, nb_clusters);
      ++ done;

      if (callback && !callback (callback_factor * (done + 1) / double(size)))
        return (nb_clusters + 1);

      neighbor_query.get_iterators (get (point_map, *current), 0, neighbor_radius,
                                    boost::make_function_output_iterator
                                    ([&](const iterator& it) { todo.push(it); }), true);

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
                                    std::back_inserter (neighbors), 0);

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
