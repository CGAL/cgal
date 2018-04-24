// Copyright (c) 2013 INRIA Sophia-Antipolis (France).
// Copyright (c) 2016 GeometryFactory Sarl (France).
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

#include <CGAL/boost/graph/named_function_params.h>
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
  struct Pmap_unary_function : public CGAL::unary_function<ValueType, typename Kernel::Point_3>
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
  struct Pmap_unary_function : public CGAL::unary_function<ValueType, typename Kernel::Point_2>
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
   \ingroup PkgPointSetProcessingAlgorithms

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
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with
     value type `geom_traits::Point_3` (or `geom_traits::Point_2`).
     If this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` (or
     `CGAL::Identity_property_map<geom_traits::Point_2>`) is
     used.\cgalParamEnd
     \cgalParamBegin{query_point_map} a model of `ReadablePropertyMap` with
     value type `geom_traits::Point_3` (or `geom_traits::Point_2`).
     If this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` (or
     `CGAL::Identity_property_map<geom_traits::Point_2>`) is
     used.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
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
  using boost::choose_param;
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::const_type PointMap;
  typedef typename Point_set_processing_3::GetQueryPointMap<QueryPointRange, NamedParameters>::const_type QueryPointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  typedef typename boost::property_traits<PointMap>::value_type Point_d;

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  QueryPointMap query_point_map = choose_param(get_param(np, internal_np::query_point_map), QueryPointMap());

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
   \ingroup PkgPointSetProcessingAlgorithms

   Estimates the global scale in a K nearest neighbors sense. The
   computed scale corresponds to the smallest scale such that the K
   subsets of points have the appearance of a surface in 3D or the
   appearance of a curve in 2D (see \ref Point_set_processing_3Scale).

   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with
     value type `geom_traits::Point_3` (or `geom_traits::Point_2`).
     If this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` (or
     `CGAL::Identity_property_map<geom_traits::Point_2>`) is
     used.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
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
  using boost::choose_param;
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::const_type PointMap;
  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
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
   \ingroup PkgPointSetProcessingAlgorithms

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
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with
     value type `geom_traits::Point_3` (or `geom_traits::Point_2`).
     If this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` (or
     `CGAL::Identity_property_map<geom_traits::Point_2>`) is
     used.\cgalParamEnd
     \cgalParamBegin{query_point_map} a model of `ReadablePropertyMap` with
     value type `geom_traits::Point_3` (or `geom_traits::Point_2`).
     If this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` (or
     `CGAL::Identity_property_map<geom_traits::Point_2>`) is
     used.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
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
  using boost::choose_param;
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::const_type PointMap;
  typedef typename Point_set_processing_3::GetQueryPointMap<QueryPointRange, NamedParameters>::const_type QueryPointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  typedef typename boost::property_traits<PointMap>::value_type Point_d;

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  QueryPointMap query_point_map = choose_param(get_param(np, internal_np::query_point_map), QueryPointMap());

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
   \ingroup PkgPointSetProcessingAlgorithms

   Estimates the global scale in a range sense. The computed scale
   corresponds to the smallest scale such that the subsets of points
   inside the sphere range have the appearance of a surface in 3D or
   the appearance of a curve in 2D (see \ref Point_set_processing_3Scale).


   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with
     value type `geom_traits::Point_3` (or `geom_traits::Point_2`).
     If this parameter is omitted,
     `CGAL::Identity_property_map<geom_traits::Point_3>` (or
     `CGAL::Identity_property_map<geom_traits::Point_2>`) is
     used.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \note This function accepts both 2D and 3D points.

   \return The estimated scale in the range sense. The return type `FT` is a number type. It is
   either deduced from the `geom_traits` \ref psp_namedparameters "Named Parameters" if provided,
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
  using boost::choose_param;
  std::vector<double> scales;
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::const_type PointMap;
  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
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

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API  
template <typename SamplesInputIterator,
          typename SamplesPointMap,
          typename QueriesInputIterator,
          typename QueriesPointMap,
          typename OutputIterator,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_local_k_neighbor_scales(), please update your code")
OutputIterator
estimate_local_k_neighbor_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  SamplesPointMap samples_map, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  QueriesPointMap queries_map, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  OutputIterator output, ///< output iterator to store the computed scales
  const Kernel& /*kernel*/) ///< geometric traits.
{
  return estimate_local_k_neighbor_scales
    (CGAL::make_range (first, beyond),
     CGAL::make_range (first_query, beyond_query),
     output,
     CGAL::parameters::point_map (samples_map).
     query_point_map (queries_map).
     geom_traits (Kernel()));
}

// deprecated API  
template <typename SamplesInputIterator,
          typename SamplesPointMap,
          typename QueriesInputIterator,
          typename QueriesPointMap,
          typename OutputIterator
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_local_k_neighbor_scales(), please update your code")
OutputIterator
estimate_local_k_neighbor_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  SamplesPointMap samples_map, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  QueriesPointMap queries_map, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  OutputIterator output) ///< output iterator to store the computed scales
{
  return estimate_local_k_neighbor_scales
    (CGAL::make_range (first, beyond),
     CGAL::make_range (first_query, beyond_query),
     output,
     CGAL::parameters::point_map (samples_map).
     query_point_map (queries_map));
}

// deprecated API
template <typename SamplesInputIterator,
          typename QueriesInputIterator,
          typename OutputIterator
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_local_k_neighbor_scales(), please update your code")
OutputIterator
estimate_local_k_neighbor_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  OutputIterator output) ///< output iterator to store the computed scales
{
  return estimate_local_k_neighbor_scales
    (CGAL::make_range (first, beyond),
     CGAL::make_range (first_query, beyond_query),
     output);
}

// deprecated API
template <typename InputIterator,
          typename PointMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_global_k_neighbor_scale(), please update your code")
std::size_t
estimate_global_k_neighbor_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  const Kernel& kernel) ///< geometric traits.
{
  return estimate_global_k_neighbor_scale
    (CGAL::make_range (first, beyond),
     CGAL::parameters::point_map (point_map).
     geom_traits (kernel));
}

// deprecated API
template <typename InputIterator,
          typename PointMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_global_k_neighbor_scale(), please update your code")
std::size_t
estimate_global_k_neighbor_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map) ///< property map: value_type of InputIterator -> Point_3 or Point_2
{
  return estimate_global_k_neighbor_scale
    (CGAL::make_range (first, beyond),
     CGAL::parameters::point_map (point_map));
}

// deprecated API  
template <typename InputIterator
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_global_k_neighbor_scale(), please update your code")
std::size_t
estimate_global_k_neighbor_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond) ///< past-the-end iterator over the input points.
{
  return estimate_global_k_neighbor_scale
    (CGAL::make_range (first, beyond));
}

// deprecated API
template <typename SamplesInputIterator,
          typename SamplesPointMap,
          typename QueriesInputIterator,
          typename QueriesPointMap,
          typename OutputIterator,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_local_range_scales(), please update your code")
OutputIterator
estimate_local_range_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  SamplesPointMap samples_map, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  QueriesPointMap queries_map, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  OutputIterator output, ///< output iterator to store the computed scales
  const Kernel& /*kernel*/) ///< geometric traits.
{
  return estimate_local_range_scales
    (CGAL::make_range (first, beyond),
     CGAL::make_range (first_query, beyond_query),
     output,
     CGAL::parameters::point_map (samples_map).
     query_point_map (queries_map).
     geom_traits (Kernel()));
}

// deprecated API
template <typename SamplesInputIterator,
          typename SamplesPointMap,
          typename QueriesInputIterator,
          typename QueriesPointMap,
          typename OutputIterator
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_local_range_scales(), please update your code")
OutputIterator
estimate_local_range_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  SamplesPointMap samples_map, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  QueriesPointMap queries_map, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  OutputIterator output) ///< output iterator to store the computed scales
{
  return estimate_local_range_scales
    (CGAL::make_range (first, beyond),
     CGAL::make_range (first_query, beyond_query),
     output,
     CGAL::parameters::point_map (samples_map).
     query_point_map (queries_map));
}

// deprecated API
template <typename SamplesInputIterator,
          typename QueriesInputIterator,
          typename OutputIterator
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_local_range_scales(), please update your code")
OutputIterator
estimate_local_range_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  OutputIterator output) ///< output iterator to store the computed scales
{
  return estimate_local_range_scales
    (CGAL::make_range (first, beyond),
     CGAL::make_range (first_query, beyond_query),
     output);
}


// deprecated API
template <typename InputIterator,
          typename PointMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_global_range_scale(), please update your code")
double
estimate_global_range_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of InputIterator -> Point_3 or Point_3
  const Kernel& kernel) ///< geometric traits.
{
  return estimate_global_range_scale
    (CGAL::make_range (first, beyond),
     CGAL::parameters::point_map (point_map).
     geom_traits (kernel));
}

// deprecated API
template <typename InputIterator,
          typename PointMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_global_range_scale(), please update your code")
double
estimate_global_range_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map) ///< property map: value_type of InputIterator -> Point_3 or Point_3
{
  return estimate_global_range_scale
    (CGAL::make_range (first, beyond),
     CGAL::parameters::point_map (point_map));
}


// deprecated API
template <typename InputIterator>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::estimate_global_range_scale(), please update your code")
double
estimate_global_range_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond) ///< past-the-end iterator over the input points.
{
  return estimate_global_range_scale
    (CGAL::make_range (first, beyond));
}
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond  

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_ESTIMATE_SCALE_3_H
