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
//
// Author(s) : Simon Giraudot

#ifndef CGAL_ESTIMATE_SCALE_H
#define CGAL_ESTIMATE_SCALE_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/assertions.h>
#include <CGAL/hierarchy_simplify_point_set.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/Point_set_2.h>

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

  template <typename ValueType, typename PointPMap>
  struct Pmap_unary_function : public std::unary_function<ValueType, typename Kernel::Point_3>
  {
    PointPMap point_pmap;
    Pmap_unary_function (PointPMap point_pmap) : point_pmap (point_pmap) { }
    const typename Kernel::Point_3& operator() (const ValueType& v) const { return get(point_pmap, v); }
  };
  
  std::size_t m_cluster_size;
  std::vector<Tree*> m_trees;
  std::vector<FT> m_weights;
  std::vector<FT> m_precomputed_factor;

public:

  template <typename InputIterator, typename PointPMap>
  Quick_multiscale_approximate_knn_distance (InputIterator first,
                                             InputIterator beyond,
                                             PointPMap point_pmap,
                                             std::size_t cluster_size = 25)
    : m_cluster_size (cluster_size)
  {
    typedef Pmap_unary_function<typename std::iterator_traits<InputIterator>::value_type,
                                PointPMap> Unary_f;

    m_trees.push_back (new Tree (boost::make_transform_iterator (first, Unary_f(point_pmap)),
                                 boost::make_transform_iterator (beyond, Unary_f(point_pmap))));
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

    InputIterator first_unused = beyond;

    nb_pts = m_trees[0]->size();
    for (std::size_t i = 1; i < nb_trees; ++ i)
      {
        first_unused
          = CGAL::hierarchy_simplify_point_set (first, first_unused, point_pmap,
                                                static_cast<unsigned int>(m_cluster_size), 1./3.);

        m_trees.push_back (new Tree(boost::make_transform_iterator (first, Unary_f(point_pmap)),
                                    boost::make_transform_iterator (first_unused, Unary_f(point_pmap))));

        m_weights.push_back (m_trees[0]->size() / (FT)(m_trees.back()->size()));
      }
  }

  ~Quick_multiscale_approximate_knn_distance()
  {
    for (std::size_t i = 0; i < m_trees.size(); ++ i)
      delete m_trees[i];
  }

  template <typename InputIterator, typename PointPMap>
  std::size_t compute_k_scale (InputIterator query, PointPMap point_pmap)
  {
    std::size_t out;
    FT dummy;
    compute_scale (query, point_pmap, out, dummy);
    return out;
  }

  template <typename InputIterator, typename PointPMap>
  FT compute_range_scale (InputIterator query, PointPMap point_pmap)
  {
    std::size_t dummy;
    FT out;
    compute_scale (query, point_pmap, dummy, out);
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
  
  
  template <typename InputIterator, typename PointPMap>
  void compute_scale (InputIterator query, PointPMap point_pmap,
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
        Neighbor_search search (*(m_trees[t]), get(point_pmap, *query),
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

  template <typename ValueType, typename PointPMap>
  struct Pmap_unary_function : public std::unary_function<ValueType, typename Kernel::Point_2>
  {
    PointPMap point_pmap;
    Pmap_unary_function (PointPMap point_pmap) : point_pmap (point_pmap) { }
    const typename Kernel::Point_2& operator() (const ValueType& v) const { return get(point_pmap, v); }
  };

  template <typename PointPMap>
  struct Pmap_to_3d
  {
    PointPMap point_pmap;
    typedef typename Kernel::Point_3 value_type;
    typedef const value_type& reference;
    typedef typename boost::property_traits<PointPMap>::key_type key_type;
    typedef boost::lvalue_property_map_tag category;
    Pmap_to_3d () { }
    Pmap_to_3d (PointPMap point_pmap)
      : point_pmap (point_pmap) { }
    friend inline value_type get (const Pmap_to_3d& ppmap, key_type i) 
    {
      typename Kernel::Point_2 p2 = get(ppmap.point_pmap, i);
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

  template <typename InputIterator, typename PointPMap>
  Quick_multiscale_approximate_knn_distance (InputIterator first,
                                             InputIterator beyond,
                                             PointPMap point_pmap,
                                             std::size_t cluster_size = 25)
    : m_cluster_size (cluster_size)
  {
    typedef Pmap_unary_function<typename std::iterator_traits<InputIterator>::value_type,
                                PointPMap> Unary_f;

    m_point_sets.push_back (new Point_set (boost::make_transform_iterator (first, Unary_f(point_pmap)),
                                           boost::make_transform_iterator (beyond, Unary_f(point_pmap))));
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

    InputIterator first_unused = beyond;
    nb_pts = m_point_sets[0]->number_of_vertices();

    for (std::size_t i = 1; i < nb_trees; ++ i)
      {
        first_unused
          = CGAL::hierarchy_simplify_point_set (first, first_unused, Pmap_to_3d<PointPMap> (point_pmap),
                                                static_cast<unsigned int>(m_cluster_size), 1./3.);

        m_point_sets.push_back (new Point_set (boost::make_transform_iterator (first, Unary_f(point_pmap)),
                                               boost::make_transform_iterator (first_unused, Unary_f(point_pmap))));

        m_weights.push_back (nb_pts / (FT)(m_point_sets.back()->number_of_vertices()));
      }

    m_cluster_size = cluster_size;
  }

  ~Quick_multiscale_approximate_knn_distance()
  {
    for (std::size_t i = 0; i < m_point_sets.size(); ++ i)
      delete m_point_sets[i];
  }

  template <typename InputIterator, typename PointPMap>
  std::size_t compute_k_scale (InputIterator query, PointPMap point_pmap)
  {
    std::size_t out;
    FT dummy;
    compute_scale (query, point_pmap, out, dummy);
    return out;
  }

  template <typename InputIterator, typename PointPMap>
  FT compute_range_scale (InputIterator query, PointPMap point_pmap)
  {
    std::size_t dummy;
    FT out;
    compute_scale (query, point_pmap, dummy, out);
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
  
  template <typename InputIterator, typename PointPMap>
  void compute_scale (InputIterator query, PointPMap point_pmap,
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
    
    const typename Kernel::Point_2& pquery = get(point_pmap, *query);
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

/// \ingroup PkgPointSetProcessingAlgorithms

/// Estimates the local scale in a K nearest neighbors sense on a set
/// of user-defined query points. The computed scales correspond to
/// the smallest scales such that the K subsets of points have the
/// appearance of a surface in 3D or the appearance of a curve in 2D
/// (see \ref Point_set_processing_3Scale).
///
///
/// @tparam SamplesInputIterator iterator over input sample points.
/// @tparam SamplesPointPMap is a model of `ReadablePropertyMap` with
///        value type `Point_3<Kernel>` or `Point_2<Kernel>`.  It can
///        be omitted if the value type of `SamplesInputIterator` is
///        convertible to `Point_3<Kernel>` or to `Point_2<Kernel>`.
/// @tparam QueriesInputIterator iterator over points where scale
///        should be computed.
/// @tparam QueriesInputIterator is a model of `ReadablePropertyMap`
///        with value type `Point_3<Kernel>` or `Point_2<Kernel>`.  It
///        can be omitted if the value type of `QueriesInputIterator` is
///        convertible to `Point_3<Kernel>` or to `Point_2<Kernel>`.
/// @tparam OutputIterator is used to store the computed scales. It accepts
///        values of type `std::size_t`.
/// @tparam Kernel Geometric traits class.  It can be omitted and
///        deduced automatically from the value type of `SamplesPointPMap`.
///
/// @note This function accepts both 2D and 3D points, but sample
///      points and query must have the same dimension.

// This variant requires all parameters.
template <typename SamplesInputIterator,
          typename SamplesPointPMap,
          typename QueriesInputIterator,
          typename QueriesPointPMap,
          typename OutputIterator,
          typename Kernel
>
OutputIterator
estimate_local_k_neighbor_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  SamplesPointPMap samples_pmap, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  QueriesPointPMap queries_pmap, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  OutputIterator output, ///< output iterator to store the computed scales
  const Kernel& /*kernel*/) ///< geometric traits.
{
  typedef typename boost::property_traits<SamplesPointPMap>::value_type Point_d;

  // Build multi-scale KD-tree
  internal::Quick_multiscale_approximate_knn_distance<Kernel, Point_d> kdtree (first, beyond, samples_pmap);

  // Compute local scales everywhere
  for (QueriesInputIterator it = first_query; it != beyond_query; ++ it)
    *(output ++) = kdtree.compute_k_scale (it, queries_pmap);

  return output;
}

  
/// \ingroup PkgPointSetProcessingAlgorithms

/// Estimates the global scale in a K nearest neighbors sense. The
/// computed scale corresponds to the smallest scale such that the K
/// subsets of points have the appearance of a surface in 3D or the
/// appearance of a curve in 2D (see \ref Point_set_processing_3Scale).
///
///
/// @tparam InputIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with
///        value type `Point_3<Kernel>` or `Point_2<Kernel>`.  It can
///        be omitted if the value type of `InputIterator` is
///        convertible to `Point_3<Kernel>` or to `Point_2<Kernel>`.
/// @tparam Kernel Geometric traits class.  It can be omitted and
///        deduced automatically from the value type of `PointPMap`.
///
/// @note This function accepts both 2D and 3D points.
///
/// @return The estimated scale in the K nearest neighbors sense.
// This variant requires all parameters.
template <typename InputIterator,
          typename PointPMap,
          typename Kernel
>
std::size_t
estimate_global_k_neighbor_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  const Kernel& kernel) ///< geometric traits.
{
  std::vector<std::size_t> scales;
  estimate_local_k_neighbor_scales (first, beyond, point_pmap,
                                    first, beyond, point_pmap,
                                    std::back_inserter (scales),
                                    kernel);
  std::sort (scales.begin(), scales.end());
  return scales[scales.size() / 2];
}

  
/// \ingroup PkgPointSetProcessingAlgorithms

/// Estimates the local scale in a range sense on a set of
/// user-defined query points. The computed scales correspond to the
/// smallest scales such that the subsets of points included in the
/// sphere range have the appearance of a surface in 3D or the
/// appearance of a curve in 2D (see \ref Point_set_processing_3Scale).
///
///
/// @tparam SamplesInputIterator iterator over input sample points.
/// @tparam SamplesPointPMap is a model of `ReadablePropertyMap` with
///        value type `Point_3<Kernel>` or `Point_2<Kernel>`.  It can
///        be omitted if the value type of `SamplesInputIterator` is
///        convertible to `Point_3<Kernel>` or to `Point_2<Kernel>`.
/// @tparam QueriesInputIterator iterator over points where scale
///        should be computed.
/// @tparam QueriesInputIterator is a model of `ReadablePropertyMap`
///        with value type `Point_3<Kernel>` or `Point_2<Kernel>`.  It
///        can be omitted if the value type of `QueriesInputIterator` is
///        convertible to `Point_3<Kernel>` or to `Point_2<Kernel>`.
/// @tparam OutputIterator is used to store the computed scales. It accepts
///        values of type `Kernel::FT`.
/// @tparam Kernel Geometric traits class.  It can be omitted and
///        deduced automatically from the value type of `SamplesPointPMap`.
///
/// @note This function accepts both 2D and 3D points, but sample
///      points and query must have the same dimension.

// This variant requires all parameters.
template <typename SamplesInputIterator,
          typename SamplesPointPMap,
          typename QueriesInputIterator,
          typename QueriesPointPMap,
          typename OutputIterator,
          typename Kernel
>
OutputIterator
estimate_local_range_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  SamplesPointPMap samples_pmap, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  QueriesPointPMap queries_pmap, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  OutputIterator output, ///< output iterator to store the computed scales
  const Kernel& /*kernel*/) ///< geometric traits.
{
  typedef typename boost::property_traits<SamplesPointPMap>::value_type Point_d;

  // Build multi-scale KD-tree
  internal::Quick_multiscale_approximate_knn_distance<Kernel, Point_d> kdtree (first, beyond, samples_pmap);

  // Compute local scales everywhere
  for (QueriesInputIterator it = first_query; it != beyond_query; ++ it)
    *(output ++) = kdtree.compute_range_scale (it, queries_pmap);

  return output;
}

  
/// \ingroup PkgPointSetProcessingAlgorithms

/// Estimates the global scale in a range sense. The computed scale
/// corresponds to the smallest scale such that the subsets of points
/// inside the sphere range have the appearance of a surface in 3D or
/// the appearance of a curve in 2D (see \ref Point_set_processing_3Scale).
///
///
/// @tparam InputIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with
///        value type `Point_3<Kernel>` or `Point_2<Kernel>`.  It can
///        be omitted if the value type of `InputIterator` is
///        convertible to `Point_3<Kernel>` or to `Point_2<Kernel>`.
/// @tparam Kernel Geometric traits class.  It can be omitted and
///        deduced automatically from the value type of `PointPMap`.
///
/// @note This function accepts both 2D and 3D points.
///
/// @return The estimated scale in the range sense.
// This variant requires all parameters.
template <typename InputIterator,
          typename PointPMap,
          typename Kernel
>
typename Kernel::FT
estimate_global_range_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3 or Point_3
  const Kernel& kernel) ///< geometric traits.
{
  std::vector<typename Kernel::FT> scales;
  estimate_local_range_scales (first, beyond, point_pmap,
                               first, beyond, point_pmap,
                               std::back_inserter (scales),
                               kernel);
  std::sort (scales.begin(), scales.end());
  return std::sqrt (scales[scales.size() / 2]);
}


// ----------------------------------------------------------------------------
// Useful overloads
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

template <typename SamplesInputIterator,
          typename SamplesPointPMap,
          typename QueriesInputIterator,
          typename QueriesPointPMap,
          typename OutputIterator
>
OutputIterator
estimate_local_k_neighbor_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  SamplesPointPMap samples_pmap, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  QueriesPointPMap queries_pmap, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  OutputIterator output) ///< output iterator to store the computed scales
{
  typedef typename boost::property_traits<SamplesPointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return estimate_local_k_neighbor_scales (first, beyond, samples_pmap, first_query, beyond_query,
                                           queries_pmap, output, Kernel());
}

template <typename SamplesInputIterator,
          typename QueriesInputIterator,
          typename OutputIterator
>
OutputIterator
estimate_local_k_neighbor_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  OutputIterator output) ///< output iterator to store the computed scales
{
  return estimate_local_k_neighbor_scales
    (first, beyond,
     make_identity_property_map (typename std::iterator_traits<SamplesInputIterator>::value_type()),
     first_query, beyond_query,
     make_identity_property_map (typename std::iterator_traits<QueriesInputIterator>::value_type()),
     output);
 }


template <typename InputIterator,
          typename PointPMap
>
std::size_t
estimate_global_k_neighbor_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap) ///< property map: value_type of InputIterator -> Point_3 or Point_2
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return estimate_global_k_neighbor_scale (first, beyond, point_pmap, Kernel());
}

template <typename InputIterator
>
std::size_t
estimate_global_k_neighbor_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond) ///< past-the-end iterator over the input points.
{
  return estimate_global_k_neighbor_scale
    (first, beyond, make_identity_property_map (typename std::iterator_traits<InputIterator>::value_type()));
}


template <typename SamplesInputIterator,
          typename SamplesPointPMap,
          typename QueriesInputIterator,
          typename QueriesPointPMap,
          typename OutputIterator
>
OutputIterator
estimate_local_range_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  SamplesPointPMap samples_pmap, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  QueriesPointPMap queries_pmap, ///< property map: value_type of InputIterator -> Point_3 or Point_2
  OutputIterator output) ///< output iterator to store the computed scales
{
  typedef typename boost::property_traits<SamplesPointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return estimate_local_range_scales(first, beyond, samples_pmap, first_query, beyond_query,
                                     queries_pmap, output, Kernel());
}


template <typename SamplesInputIterator,
          typename QueriesInputIterator,
          typename OutputIterator
>
OutputIterator
estimate_local_range_scales(
  SamplesInputIterator first, ///< iterator over the first input sample.
  SamplesInputIterator beyond, ///< past-the-end iterator over the input samples.
  QueriesInputIterator first_query, ///< iterator over the first point where scale must be estimated
  QueriesInputIterator beyond_query, ///< past-the-end iterator over the points where scale must be estimated
  OutputIterator output) ///< output iterator to store the computed scales
{
  return estimate_local_range_scales
    (first, beyond,
     make_identity_property_map (typename std::iterator_traits<SamplesInputIterator>::value_type()),
     first_query, beyond_query,
     make_identity_property_map (typename std::iterator_traits<QueriesInputIterator>::value_type()),
     output);
}



template <typename InputIterator,
          typename PointPMap
>
typename Kernel_traits<typename boost::property_traits<PointPMap>::value_type>::Kernel::FT
estimate_global_range_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap) ///< property map: value_type of InputIterator -> Point_3 or Point_3
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return estimate_global_range_scale (first, beyond, point_pmap, Kernel());
}



template <typename InputIterator>
typename Kernel_traits<typename std::iterator_traits<InputIterator>::value_type>::Kernel::FT
estimate_global_range_scale(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond) ///< past-the-end iterator over the input points.
{
  return estimate_global_range_scale
    (first, beyond, make_identity_property_map (typename std::iterator_traits<InputIterator>::value_type()));
                                      
}
/// \endcond  

} //namespace CGAL

#endif // CGAL_ESTIMATE_SCALE_3_H
