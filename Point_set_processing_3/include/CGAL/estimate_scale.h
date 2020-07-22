// Copyright (c) 2013 INRIA Sophia-Antipolis (France).
// Copyright (c) 2016 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_ESTIMATE_SCALE_H
#define CGAL_ESTIMATE_SCALE_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/assertions.h>
#include <CGAL/hierarchy_simplify_point_set.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/Point_set_2.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <fstream>

#include <iterator>
#include <list>


namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {

template <class Kernel, class PointType>
class Quick_multiscale_approximate_knn_distance
{

};


template <class Kernel>
class Quick_multiscale_approximate_knn_distance<Kernel, typename Kernel::Point_3>
{
  typedef typename Kernel::FT FT;
  typedef Search_traits_3<Kernel> Tree_traits;
  typedef Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Iterator;

  template <typename ValueType, typename PointMap>
  struct Pmap_unary_function : public CGAL::cpp98::unary_function<ValueType, typename Kernel::Point_3>
  {
    PointMap point_map;
    Pmap_unary_function (PointMap point_map) : point_map (point_map) { }
    typename boost::property_traits<PointMap>::reference
    operator() (const ValueType& v) const { return get(point_map, v); }
  };

  std::size_t m_cluster_size;
  std::vector<Tree*> m_trees;
  std::vector<FT> m_weights;
  std::vector<FT> m_precomputed_factor;

public:

  template <typename InputIterator, typename PointMap>
  Quick_multiscale_approximate_knn_distance (InputIterator first,
                                             InputIterator beyond,
                                             PointMap point_map,
                                             std::size_t cluster_size = 25)
    : m_cluster_size (cluster_size)
  {
    typedef Pmap_unary_function<typename std::iterator_traits<InputIterator>::value_type,
                                PointMap> Unary_f;

    // Avoid moving points of input as the range is const
    std::vector<typename Kernel::Point_3> kd_tree_points;
    std::copy (boost::make_transform_iterator (first, Unary_f(point_map)),
               boost::make_transform_iterator (beyond, Unary_f(point_map)),
               std::back_inserter (kd_tree_points));


    m_trees.push_back (new Tree (kd_tree_points.begin(), kd_tree_points.end()));
    m_weights.push_back (1.);
    std::size_t nb_pts = m_trees[0]->size();

    std::size_t nb_trees = 0;
    while (nb_pts > m_cluster_size)
      {
        nb_trees ++;
        nb_pts /= m_cluster_size;
      }

    m_trees.reserve (nb_trees);
    m_weights.reserve (nb_trees);

    typename std::vector<typename Kernel::Point_3>::iterator first_unused = kd_tree_points.end();

    nb_pts = m_trees[0]->size();
    for (std::size_t i = 1; i < nb_trees; ++ i)
      {
        CGAL::Iterator_range<typename std::vector<typename Kernel::Point_3>::iterator> points
          (kd_tree_points.begin(), first_unused);
        first_unused
          = CGAL::hierarchy_simplify_point_set (points,
                                                CGAL::parameters::size(static_cast<unsigned int>(m_cluster_size)).
                                                maximum_variation(1./3.));

        m_trees.push_back (new Tree(kd_tree_points.begin(), first_unused));

        m_weights.push_back (m_trees[0]->size() / (FT)(m_trees.back()->size()));
      }
  }

  ~Quick_multiscale_approximate_knn_distance()
  {
    for (std::size_t i = 0; i < m_trees.size(); ++ i)
      delete m_trees[i];
  }

  template <typename InputIterator, typename PointMap>
  std::size_t compute_k_scale (InputIterator query, PointMap point_map)
  {
    std::size_t out;
    FT dummy;
    compute_scale (query, point_map, out, dummy);
    return out;
  }

  template <typename InputIterator, typename PointMap>
  FT compute_range_scale (InputIterator query, PointMap point_map)
  {
    std::size_t dummy;
    FT out;
    compute_scale (query, point_map, dummy, out);
    return out;
  }

  void precompute_factors ()
  {
    FT nb = 0.;
    for (std::size_t t = 0; t < m_trees.size(); ++ t)
      {
        std::size_t size = (t == (m_trees.size() - 1)
                            ? m_trees[t]->size()
                            : static_cast<std::size_t>(m_weights[t+1] / m_weights[t]));
        for (std::size_t i = (t == 0 ? 0 : 1); i < size; ++ i)
          {
            nb += m_weights[t];
            if (nb < 6.) // do not consider values under 6
              continue;
            m_precomputed_factor.push_back (0.91666666 * std::log (nb));
          }
      }
  }


  template <typename InputIterator, typename PointMap>
  void compute_scale (InputIterator query, PointMap point_map,
                      std::size_t& k, FT& d)
  {
    if (m_precomputed_factor.empty())
      precompute_factors();

    k = 0;
    d = 0.;

    FT dist_min = (std::numeric_limits<FT>::max)();
    FT sum_sq_distances = 0.;
    FT nb = 0.;
    std::size_t index = 0;
    for (std::size_t t = 0; t < m_trees.size(); ++ t)
      {
        Neighbor_search search (*(m_trees[t]), get(point_map, *query),
                                static_cast<unsigned int>((t == (m_trees.size() - 1)
                                                           ? m_trees[t]->size()
                                                           : m_weights[t+1] / m_weights[t])));
        Iterator it = search.begin();

        if (t != 0) // Skip first point except on first scale
          ++ it;

        for (; it != search.end(); ++ it)
          {
            sum_sq_distances += m_weights[t] * it->second;
            nb += m_weights[t];

            if (nb < 6.) // do not consider values under 6
              continue;

            // sqrt(sum_sq_distances / nb) / nb^(5/12)
            // Computed in log space with precomputed factor for time optimization
            FT dist = 0.5 * std::log (sum_sq_distances) - m_precomputed_factor[index ++];

            if (dist < dist_min)
              {
                dist_min = dist;
                k = (std::size_t)nb;
                d = it->second;
              }
          }
      }
  }

};


template <class Kernel>
class Quick_multiscale_approximate_knn_distance<Kernel, typename Kernel::Point_2>
{
  typedef typename Kernel::FT FT;
  typedef CGAL::Point_set_2<Kernel> Point_set;
  typedef typename Point_set::Vertex_handle Vertex_handle;

  template <typename ValueType, typename PointMap>
  struct Pmap_unary_function : public CGAL::cpp98::unary_function<ValueType, typename Kernel::Point_2>
  {
    PointMap point_map;
    Pmap_unary_function (PointMap point_map) : point_map (point_map) { }
    typename boost::property_traits<PointMap>::reference
    operator() (const ValueType& v) const { return get(point_map, v); }
  };

  template <typename PointMap>
  struct Pmap_to_3d
  {
    PointMap point_map;
    typedef typename Kernel::Point_3 value_type;
    typedef const value_type& reference;
    typedef typename Kernel::Point_2 key_type;
    typedef boost::lvalue_property_map_tag category;

    Pmap_to_3d () { }
    Pmap_to_3d (PointMap point_map)
      : point_map (point_map) { }

    friend inline value_type get (const Pmap_to_3d& pmap, key_type p)
    {
      typename boost::property_traits<PointMap>::reference
        p2 = get(pmap.point_map, p);
      return value_type (p2.x(), p2.y(), 0.);
    }

  };

  struct Sort_by_distance_to_point
  {
    const typename Kernel::Point_2& ref;
    Sort_by_distance_to_point (const typename Kernel::Point_2& ref) : ref (ref) { }
    bool operator() (const Vertex_handle& a, const Vertex_handle& b)
    {
      return (CGAL::squared_distance (a->point(), ref)
              < CGAL::squared_distance (b->point(), ref));
    }
  };


  std::size_t m_cluster_size;
  std::vector<Point_set*> m_point_sets;
  std::vector<FT> m_weights;
  std::vector<FT> m_precomputed_factor;

public:

  template <typename InputIterator, typename PointMap>
  Quick_multiscale_approximate_knn_distance (InputIterator first,
                                             InputIterator beyond,
                                             PointMap point_map,
                                             std::size_t cluster_size = 25)
    : m_cluster_size (cluster_size)
  {
    typedef Pmap_unary_function<typename std::iterator_traits<InputIterator>::value_type,
                                PointMap> Unary_f;

    // Avoid moving points of input as the range is const
    std::vector<typename Kernel::Point_2> search_points;
    std::copy (boost::make_transform_iterator (first, Unary_f(point_map)),
               boost::make_transform_iterator (beyond, Unary_f(point_map)),
               std::back_inserter (search_points));


    m_point_sets.push_back (new Point_set (search_points.begin(), search_points.end()));
    m_weights.push_back (1.);

    std::size_t nb_pts = m_point_sets[0]->number_of_vertices();
    std::size_t nb_trees = 0;
    while (nb_pts > m_cluster_size)
      {
        nb_trees ++;
        nb_pts /= m_cluster_size;
      }

    m_point_sets.reserve (nb_trees);
    m_weights.reserve (nb_trees);

    typename std::vector<typename Kernel::Point_2>::iterator first_unused = search_points.end();
    nb_pts = m_point_sets[0]->number_of_vertices();

    for (std::size_t i = 1; i < nb_trees; ++ i)
      {
        CGAL::Iterator_range<typename std::vector<typename Kernel::Point_2>::iterator> points
          (search_points.begin(), first_unused);
        first_unused
          = CGAL::hierarchy_simplify_point_set (points,
                                                CGAL::parameters::point_map(Pmap_to_3d<PointMap>(point_map)).
                                                size(static_cast<unsigned int>(m_cluster_size)).
                                                maximum_variation(1./3.));

        m_point_sets.push_back (new Point_set (search_points.begin(), first_unused));

        m_weights.push_back (nb_pts / (FT)(m_point_sets.back()->number_of_vertices()));
      }

    m_cluster_size = cluster_size;
  }

  ~Quick_multiscale_approximate_knn_distance()
  {
    for (std::size_t i = 0; i < m_point_sets.size(); ++ i)
      delete m_point_sets[i];
  }

  template <typename InputIterator, typename PointMap>
  std::size_t compute_k_scale (InputIterator query, PointMap point_map)
  {
    std::size_t out;
    FT dummy;
    compute_scale (query, point_map, out, dummy);
    return out;
  }

  template <typename InputIterator, typename PointMap>
  FT compute_range_scale (InputIterator query, PointMap point_map)
  {
    std::size_t dummy;
    FT out;
    compute_scale (query, point_map, dummy, out);
    return out;
  }

  void precompute_factors ()
  {
    FT nb = 0.;
    for (std::size_t t = 0; t < m_point_sets.size(); ++ t)
      {
        std::size_t size = (t == m_point_sets.size() - 1
                            ? m_point_sets[t]->number_of_vertices()
                            : static_cast<std::size_t>(m_weights[t+1] / m_weights[t]));
        for (std::size_t i = (t == 0 ? 0 : 1); i < size; ++ i)
          {
            nb += m_weights[t];
            if (nb < 6.) // do not consider values under 6
              continue;
            m_precomputed_factor.push_back (1.25 * std::log (nb));
          }
      }
  }

  template <typename InputIterator, typename PointMap>
  void compute_scale (InputIterator query, PointMap point_map,
                      std::size_t& k, FT& d)
  {
    if (m_precomputed_factor.empty())
      precompute_factors();

    k = 0;
    d = 0.;

    FT dist_min = (std::numeric_limits<FT>::max)();
    FT sum_sq_distances = 0.;
    FT nb = 0.;
    std::size_t index = 0;

    typename boost::property_traits<PointMap>::reference
      pquery = get(point_map, *query);
    for (std::size_t t = 0; t < m_point_sets.size(); ++ t)
      {
        std::size_t size = ((t == m_point_sets.size() - 1)
                            ? m_point_sets[t]->number_of_vertices()
                            : static_cast<std::size_t>(m_weights[t+1] / m_weights[t]));
        std::vector<Vertex_handle> neighbors;
        neighbors.reserve (size);
        m_point_sets[t]->nearest_neighbors (pquery, size, std::back_inserter (neighbors));

        std::sort (neighbors.begin(), neighbors.end(),
                   Sort_by_distance_to_point (pquery));
        for (std::size_t n = (t == 0 ? 0 : 1); n < neighbors.size(); ++ n)
          {
            FT sq_dist = CGAL::squared_distance (pquery, neighbors[n]->point());

            sum_sq_distances += m_weights[t] * sq_dist;
            nb += m_weights[t];

            if (nb < 6.) // do not consider values under 6
              continue;

            // sqrt(sum_sq_distances / nb) / nb^(3/4)
            // Computed in log space with precomputed factor for time optimization
            FT dist = 0.5 * std::log (sum_sq_distances) - m_precomputed_factor[index ++];

            if (dist < dist_min)
              {
                dist_min = dist;
                k = (std::size_t)nb;
                d = sq_dist;
              }
          }
      }
  }

};

} /* namespace internal */
/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessing3Algorithms

   Estimates the local scale in a K nearest neighbors sense on a set
   of user-defined query points. The computed scales correspond to
   the smallest scales such that the K subsets of points have the
   appearance of a surface in 3D or the appearance of a curve in 2D
   (see \ref Point_set_processing_3Scale).

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam QueryPointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `query_point_map`.
   \tparam OutputIterator is used to store the computed scales. It accepts
   values of type `std::size_t`.

   \param points input point range.
   \param queries range of locations where scale must be estimated
   \param output iterator to store the computed scales
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`
                      (or `geom_traits::Point_2`)}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>` (or
                         `CGAL::Identity_property_map<geom_traits::Point_2>`)}
     \cgalParamNEnd

     \cgalParamNBegin{query_point_map}
       \cgalParamDescription{the property map containing the points associated to the elements of the point range `queries`}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3` (or `geom_traits::Point_2`)}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>` (or
                         `CGAL::Identity_property_map<geom_traits::Point_2>`)}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \note This function accepts both 2D and 3D points, but sample
   points and query must have the same dimension.
*/
template <typename PointRange,
          typename QueryPointRange,
          typename OutputIterator,
          typename NamedParameters
>
OutputIterator
estimate_local_k_neighbor_scales(
  const PointRange& points,
  const QueryPointRange& queries,
  OutputIterator output,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::const_type PointMap;
  typedef typename Point_set_processing_3::GetQueryPointMap<QueryPointRange, NamedParameters>::const_type QueryPointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  typedef typename boost::property_traits<PointMap>::value_type Point_d;

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  QueryPointMap query_point_map = choose_parameter<QueryPointMap>(get_parameter(np, internal_np::query_point_map));

  // Build multi-scale KD-tree
  internal::Quick_multiscale_approximate_knn_distance<Kernel, Point_d> kdtree (points.begin(),
                                                                               points.end(),
                                                                               point_map);

  // Compute local scales everywhere
  for (typename QueryPointRange::const_iterator it = queries.begin();
       it != queries.end(); ++ it)
    *(output ++) = kdtree.compute_k_scale (it, query_point_map);

  return output;
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange,
          typename QueryPointRange,
          typename OutputIterator
>
OutputIterator
estimate_local_k_neighbor_scales(
  const PointRange& points,
  const QueryPointRange& queries,
  OutputIterator output)
{
  return estimate_local_k_neighbor_scales
    (points, queries, output, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond

/**
   \ingroup PkgPointSetProcessing3Algorithms

   Estimates the global scale in a K nearest neighbors sense. The
   computed scale corresponds to the smallest scale such that the K
   subsets of points have the appearance of a surface in 3D or the
   appearance of a curve in 2D (see \ref Point_set_processing_3Scale).

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`
                      (or `geom_traits::Point_2`)}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>` (or
                         `CGAL::Identity_property_map<geom_traits::Point_2>`)}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \note This function accepts both 2D and 3D points.

   \return The estimated scale in the K nearest neighbors sense.
*/
template <typename PointRange,
          typename NamedParameters
>
std::size_t
estimate_global_k_neighbor_scale(
  const PointRange& points,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::const_type PointMap;
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  std::vector<std::size_t> scales;
  estimate_local_k_neighbor_scales (points, points, std::back_inserter (scales), np.query_point_map(point_map));
  std::sort (scales.begin(), scales.end());
  return scales[scales.size() / 2];
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange>
std::size_t
estimate_global_k_neighbor_scale(const PointRange& points)
{
  return estimate_global_k_neighbor_scale
    (points, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond

/**
   \ingroup PkgPointSetProcessing3Algorithms

   Estimates the local scale in a range sense on a set of
   user-defined query points. The computed scales correspond to the
   smallest scales such that the subsets of points included in the
   sphere range have the appearance of a surface in 3D or the
   appearance of a curve in 2D (see \ref Point_set_processing_3Scale).


   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.
   \tparam QueryPointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `query_point_map`.
   \tparam OutputIterator is used to store the computed scales. It accepts
   values of type `geom_traits::FT`.

   \param points input point range.
   \param queries range of locations where scale must be estimated
   \param output iterator to store the computed scales
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`
                      (or `geom_traits::Point_2`)}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>` (or
                         `CGAL::Identity_property_map<geom_traits::Point_2>`)}
     \cgalParamNEnd

     \cgalParamNBegin{query_point_map}
       \cgalParamDescription{the property map containing the points associated to the elements of the point range `queries`}
       \cgalParamType{a model of `ReadablePropertyMap` with value type `geom_traits::Point_3` (or `geom_traits::Point_2`)}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>` (or
                         `CGAL::Identity_property_map<geom_traits::Point_2>`)}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \note This function accepts both 2D and 3D points, but sample
   points and query must have the same dimension.
*/
template <typename PointRange,
          typename QueryPointRange,
          typename OutputIterator,
          typename NamedParameters
>
OutputIterator
estimate_local_range_scales(
  const PointRange& points,
  const QueryPointRange& queries,
  OutputIterator output,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::const_type PointMap;
  typedef typename Point_set_processing_3::GetQueryPointMap<QueryPointRange, NamedParameters>::const_type QueryPointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  typedef typename boost::property_traits<PointMap>::value_type Point_d;

  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  QueryPointMap query_point_map = choose_parameter<QueryPointMap>(get_parameter(np, internal_np::query_point_map));

  // Build multi-scale KD-tree
  internal::Quick_multiscale_approximate_knn_distance<Kernel, Point_d> kdtree (points.begin(),
                                                                               points.end(),
                                                                               point_map);

  // Compute local scales everywhere
  for (typename QueryPointRange::const_iterator it = queries.begin(); it != queries.end(); ++ it)
    *(output ++) = kdtree.compute_range_scale (it, query_point_map);

  return output;
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange,
          typename QueryPointRange,
          typename OutputIterator
>
OutputIterator
estimate_local_range_scales(
  const PointRange& points,
  const QueryPointRange& queries,
  OutputIterator output)
{
  return estimate_local_range_scales
    (points, queries, output, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond

/**
   \ingroup PkgPointSetProcessing3Algorithms

   Estimates the global scale in a range sense. The computed scale
   corresponds to the smallest scale such that the subsets of points
   inside the sphere range have the appearance of a surface in 3D or
   the appearance of a curve in 2D (see \ref Point_set_processing_3Scale).


   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

   \cgalNamedParamsBegin
     \cgalParamNBegin{point_map}
       \cgalParamDescription{a property map associating points to the elements of the point set `points`}
       \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
                      of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`
                      (or `geom_traits::Point_2`)}
       \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>` (or
                         `CGAL::Identity_property_map<geom_traits::Point_2>`)}
     \cgalParamNEnd

     \cgalParamNBegin{geom_traits}
       \cgalParamDescription{an instance of a geometric traits class}
       \cgalParamType{a model of `Kernel`}
       \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
     \cgalParamNEnd
   \cgalNamedParamsEnd

   \note This function accepts both 2D and 3D points.

   \return The estimated scale in the range sense. The return type `FT` is a number type. It is
   either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
   or the geometric traits class deduced from the point property map
   of `points`.
*/
template <typename PointRange,
          typename NamedParameters
>
#ifdef DOXYGEN_RUNNING
  FT
#else
  typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel::FT
#endif
estimate_global_range_scale(
  const PointRange& points,
  const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  std::vector<double> scales;
  typedef typename CGAL::GetPointMap<PointRange, NamedParameters>::const_type PointMap;
  PointMap point_map = choose_parameter<PointMap>(get_parameter(np, internal_np::point_map));
  estimate_local_range_scales (points, points, std::back_inserter (scales), np.query_point_map(point_map));
  std::sort (scales.begin(), scales.end());
  return std::sqrt (scales[scales.size() / 2]);
}

/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename PointRange>
typename Point_set_processing_3::GetFT<PointRange>::type
estimate_global_range_scale(const PointRange& points)
{
  return estimate_global_range_scale
    (points, CGAL::Point_set_processing_3::parameters::all_default(points));
}
/// \endcond

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_ESTIMATE_SCALE_3_H
