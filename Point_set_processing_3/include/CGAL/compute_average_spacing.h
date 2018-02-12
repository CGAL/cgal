// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
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
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_AVERAGE_SPACING_3_H
#define CGAL_AVERAGE_SPACING_3_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/assertions.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iterator>
#include <list>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>  
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {


/// Computes average spacing of one query point from K nearest neighbors.
///
/// \pre `k >= 2`.
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return average spacing (scalar).
template < typename Kernel,
           typename Tree >
typename Kernel::FT
compute_average_spacing(const typename Kernel::Point_3& query, ///< 3D point whose spacing we want to compute
                        const Tree& tree,                      ///< KD-tree
                        unsigned int k)                        ///< number of neighbors
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;

  // types for K nearest neighbors search
  typedef Search_traits_3<Kernel> Tree_traits;
  typedef Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // performs k + 1 queries (if unique the query point is
  // output first). search may be aborted when k is greater
  // than number of input points
  Neighbor_search search(tree,query,k+1);
  Search_iterator search_iterator = search.begin();
  FT sum_distances = (FT)0.0;
  unsigned int i;
  for(i=0;i<(k+1);i++)
  {
    if(search_iterator == search.end())
      break; // premature ending

    Point p = search_iterator->first;
    sum_distances += std::sqrt(CGAL::squared_distance(query,p));
    search_iterator++;
  }

  // output average spacing
  return sum_distances / (FT)i;
}


#ifdef CGAL_LINKED_WITH_TBB
  template <typename Kernel, typename Tree>
  class Compute_average_spacings {
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::FT FT;
    const Tree& tree;
    const unsigned int k;
    const std::vector<Point>& input;
    std::vector<FT>& output;

  public:
    Compute_average_spacings(Tree& tree, unsigned int k, std::vector<Point>& points,
			     std::vector<FT>& output)
      : tree(tree), k (k), input (points), output (output)
    { }
    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for( std::size_t i = r.begin(); i != r.end(); ++i)
	output[i] = CGAL::internal::compute_average_spacing<Kernel,Tree>(input[i], tree, k);
    }

  };
#endif // CGAL_LINKED_WITH_TBB

} /* namespace internal */
/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/**
   \ingroup PkgPointSetProcessingAlgorithms
   Computes average spacing from k nearest neighbors.

   \pre `k >= 2.`

   \tparam ConcurrencyTag enables sequential versus parallel algorithm.
   Possible values are `Sequential_tag`
   and `Parallel_tag`.
   \tparam PointRange is a model of `ConstRange`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param k number of neighbors.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \return average spacing (scalar). The return type `FT` is a number type. It is
   either deduced from the `geom_traits` \ref psp_namedparameters "Named Parameters" if provided,
   or the geometric traits class deduced from the point property map
   of `points`.
*/
template <typename ConcurrencyTag,
	  typename PointRange,
          typename NamedParameters
>
#ifdef DOXYGEN_RUNNING
  FT
#else
  typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel::FT
#endif
compute_average_spacing(
  const PointRange& points,
  unsigned int k,
  const NamedParameters& np)
{
  using boost::choose_param;

  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::const_type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;

  typedef typename Kernel::Point_3 Point;

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  
  // types for K nearest neighbors search structure
  typedef typename Kernel::FT FT;
  typedef Search_traits_3<Kernel> Tree_traits;
  typedef Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  // Instanciate a KD-tree search.
  // Note: We have to convert each input iterator to Point_3.
  std::vector<Point> kd_tree_points; 
  for(typename PointRange::const_iterator it = points.begin(); it != points.end(); it++)
    kd_tree_points.push_back(get(point_map, *it));
  Tree tree(kd_tree_points.begin(), kd_tree_points.end());

  // iterate over input points, compute and output normal
  // vectors (already normalized)
  FT sum_spacings = (FT)0.0;

#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
			     "Parallel_tag is enabled but TBB is unavailable.");
#else
   if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
   {
     std::vector<FT> spacings (kd_tree_points.size ());
     CGAL::internal::Compute_average_spacings<Kernel, Tree>
       f (tree, k, kd_tree_points, spacings);
     tbb::parallel_for(tbb::blocked_range<size_t>(0, kd_tree_points.size ()), f);
     for (unsigned int i = 0; i < spacings.size (); ++ i)
       sum_spacings += spacings[i];
   }
   else
#endif
     {
       for(typename PointRange::const_iterator it = points.begin(); it != points.end(); it++)
	 {
	   sum_spacings += internal::compute_average_spacing<Kernel,Tree>(
									  get(point_map,*it),
									  tree,k);
	 }
     }
   
  // return average spacing
   return sum_spacings / (FT)(kd_tree_points.size ());
}

/// \cond SKIP_IN_MANUAL

// variant with default NP
template <typename ConcurrencyTag, typename PointRange>
typename Point_set_processing_3::GetFT<PointRange>::type
compute_average_spacing(
  const PointRange& points,
  unsigned int k) ///< number of neighbors.
{
  return compute_average_spacing<ConcurrencyTag>
    (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}


#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template <typename ConcurrencyTag,
	  typename InputIterator,
          typename PointMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::compute_average_spacing(), please update your code")
typename Kernel::FT
compute_average_spacing(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of InputIterator -> Point_3
  unsigned int k, ///< number of neighbors.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  return compute_average_spacing<ConcurrencyTag>(
    CGAL::make_range (first,beyond),
    k,
    CGAL::parameters::point_map (point_map).geom_traits (Kernel()));
}

  

// deprecated API
template <typename ConcurrencyTag,
	  typename InputIterator,
          typename PointMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::compute_average_spacing(), please update your code")
typename Kernel_traits<typename boost::property_traits<PointMap>::value_type>::Kernel::FT
compute_average_spacing(
  InputIterator first,    ///< iterator over the first input point.
  InputIterator beyond,   ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of InputIterator -> Point_3
  unsigned int k) ///< number of neighbors
{
  return compute_average_spacing<ConcurrencyTag>(
    CGAL::make_range (first,beyond),
    k,
    CGAL::parameters::point_map (point_map));
}

// deprecated API
template < typename ConcurrencyTag, typename InputIterator >
typename Kernel_traits<typename std::iterator_traits<InputIterator>::value_type>::Kernel::FT
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::compute_average_spacing(), please update your code")
compute_average_spacing(
  InputIterator first,    ///< iterator over the first input point.
  InputIterator beyond,   ///< past-the-end iterator over the input points.
  unsigned int k) ///< number of neighbors.
{
  return compute_average_spacing<ConcurrencyTag>(
    CGAL::make_range (first,beyond), k);
}
#endif // CGAL_NO_DEPRECATED_CODE
  
/// \endcond


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_AVERAGE_SPACING_3_H
