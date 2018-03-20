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
// Author(s) : Pierre Alliez, Marc Pouget and Laurent Saboret

#ifndef CGAL_JET_SMOOTH_POINT_SET_H
#define CGAL_JET_SMOOTH_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/trace.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

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

/// Smoothes one point position using jet fitting on the k
/// nearest neighbors and reprojection onto the jet.
///
/// \pre `k >= 2`
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return computed point
template <typename Kernel,
          typename SvdTraits,
          typename Tree
          >
typename Kernel::Point_3
jet_smooth_point(
  const typename Kernel::Point_3& query, ///< 3D point to project
  Tree& tree, ///< KD-tree
  const unsigned int k, ///< number of neighbors.
  const unsigned int degree_fitting,
  const unsigned int degree_monge)
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;

  // types for K nearest neighbors search
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // types for jet fitting
  typedef Monge_via_jet_fitting< Kernel,
                                 Simple_cartesian<double>,
                                 SvdTraits> Monge_jet_fitting;
  typedef typename Monge_jet_fitting::Monge_form Monge_form;

  // Gather set of (k+1) neighboring points.
  // Performs k + 1 queries (if unique the query point is
  // output first). Search may be aborted if k is greater
  // than number of input points.
  std::vector<Point> points; points.reserve(k+1);
  Neighbor_search search(tree,query,k+1);
  Search_iterator search_iterator = search.begin();
  unsigned int i;
  for(i=0;i<(k+1);i++)
  {
    if(search_iterator == search.end())
      break; // premature ending
    points.push_back(search_iterator->first);
    search_iterator++;
  }
  CGAL_point_set_processing_precondition(points.size() >= 1);

  // performs jet fitting
  Monge_jet_fitting monge_fit;
  Monge_form monge_form = monge_fit(points.begin(), points.end(),
                                    degree_fitting, degree_monge);

  // output projection of query point onto the jet
  return monge_form.origin();
}

#ifdef CGAL_LINKED_WITH_TBB
  template <typename Kernel, typename SvdTraits, typename Tree>
  class Jet_smooth_pwns {
    typedef typename Kernel::Point_3 Point;
    const Tree& tree;
    const unsigned int k;
    unsigned int degree_fitting;
    unsigned int degree_monge;
    const std::vector<Point>& input;
    std::vector<Point>& output;

  public:
    Jet_smooth_pwns (Tree& tree, unsigned int k, std::vector<Point>& points,
		     unsigned int degree_fitting, unsigned int degree_monge, std::vector<Point>& output)
      : tree(tree), k (k), degree_fitting (degree_fitting),
	degree_monge (degree_monge), input (points), output (output)
    { }
    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for( std::size_t i = r.begin(); i != r.end(); ++i)
	output[i] = CGAL::internal::jet_smooth_point<Kernel, SvdTraits>(input[i], tree, k,
									degree_fitting,
									degree_monge);
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
   Smoothes the range of `points` using jet fitting on the k
   nearest neighbors and reprojection onto the jet.
   As this method relocates the points, it
   should not be called on containers sorted w.r.t. point locations.

   \pre `k >= 2`

   \tparam ConcurrencyTag enables sequential versus parallel algorithm.
   Possible values are `Sequential_tag`
   and `Parallel_tag`.
   \tparam PointRange is a model of `Range`. The value type of
   its iterator is the key type of the named parameter `point_map`.

   \param points input point range.
   \param k number of neighbors
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `ReadablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{degree_fitting} degree of jet fitting.\cgalParamEnd
     \cgalParamBegin{degree_monge} Monge degree.\cgalParamEnd
     \cgalParamBegin{svd_traits} template parameter for the class `Monge_via_jet_fitting`. If
     \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined,
     then `CGAL::Eigen_svd` is used.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

*/
template <typename ConcurrencyTag,
	  typename PointRange,
          typename NamedParameters
>
void
jet_smooth_point_set(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np)
{
  using boost::choose_param;
  
  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
  typedef typename GetSvdTraits<NamedParameters>::type SvdTraits;

  CGAL_static_assertion_msg(!(boost::is_same<SvdTraits,
                              typename GetSvdTraits<NamedParameters>::NoTraits>::value),
                            "Error: no SVD traits");

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  unsigned int degree_fitting = choose_param(get_param(np, internal_np::degree_fitting), 2);
  unsigned int degree_monge = choose_param(get_param(np, internal_np::degree_monge), 2);

  typedef typename Kernel::Point_3 Point;

  // types for K nearest neighbors search structure
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(points.begin() != points.end());

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);
  
  typename PointRange::iterator it;

  // Instanciate a KD-tree search.
  // Note: We have to convert each input iterator to Point_3.
  std::vector<Point> kd_tree_points; 
  for(it = points.begin(); it != points.end(); it++)
    kd_tree_points.push_back(get(point_map, *it));
  Tree tree(kd_tree_points.begin(), kd_tree_points.end());

  // Iterates over input points and mutates them.
  // Implementation note: the cast to Point& allows to modify only the point's position.

#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
			     "Parallel_tag is enabled but TBB is unavailable.");
#else
   if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
   {
     std::vector<Point> mutated_points (kd_tree_points.size ());
     CGAL::internal::Jet_smooth_pwns<Kernel, SvdTraits, Tree>
       f (tree, k, kd_tree_points, degree_fitting, degree_monge,
	  mutated_points);
     tbb::parallel_for(tbb::blocked_range<size_t>(0, kd_tree_points.size ()), f);
     unsigned int i = 0;
     for(it = points.begin(); it != points.end(); ++ it, ++ i)
       {
	 put(point_map, *it, mutated_points[i]);

       }
   }
   else
#endif
     {
       for(it = points.begin(); it != points.end(); it++)
	 {
	   const typename boost::property_traits<PointMap>::reference p = get(point_map, *it);
	   put(point_map, *it ,
	       internal::jet_smooth_point<Kernel, SvdTraits>(
							     p,tree,k,degree_fitting,degree_monge) );
	 }
     }
}


/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename ConcurrencyTag,
	  typename PointRange>
void
jet_smooth_point_set(
  PointRange& points,
  unsigned int k) ///< number of neighbors.
{
  jet_smooth_point_set<ConcurrencyTag>
    (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}


#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template <typename ConcurrencyTag,
	  typename InputIterator,
          typename PointMap,
          typename Kernel,
          typename SvdTraits
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::jet_smooth_point_set(), please update your code")
void
jet_smooth_point_set(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of InputIterator -> Point_3.
  unsigned int k, ///< number of neighbors.
  const Kernel& /*kernel*/, ///< geometric traits.
  unsigned int degree_fitting = 2, ///< fitting degree
  unsigned int degree_monge = 2)  ///< Monge degree
{
  CGAL::Iterator_range<InputIterator> points (first, beyond);
  return jet_smooth_point_set<ConcurrencyTag>
    (points,
     k,
     CGAL::parameters::point_map (point_map).
     degree_fitting (degree_fitting).
     degree_monge (degree_monge).
     geom_traits(Kernel()));
}
  
#if defined(CGAL_EIGEN3_ENABLED) || defined(CGAL_LAPACK_ENABLED)
// deprecated API
template <typename ConcurrencyTag,
	  typename InputIterator,
          typename PointMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::jet_smooth_point_set(), please update your code")
void
jet_smooth_point_set(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of InputIterator -> Point_3.
  unsigned int k, ///< number of neighbors.
  const Kernel& kernel, ///< geometric traits.
  unsigned int degree_fitting = 2, ///< fitting degree
  unsigned int degree_monge = 2)  ///< Monge degree
{
  #ifdef CGAL_EIGEN3_ENABLED
  typedef Eigen_svd SvdTraits;
  #else
  typedef Lapack_svd SvdTraits;
  #endif
  CGAL::Iterator_range<InputIterator> points (first, beyond);
  return jet_smooth_point_set<ConcurrencyTag>
    (points,
     k,
     CGAL::parameters::point_map (point_map).
     degree_fitting (degree_fitting).
     degree_monge (degree_monge).
     svd_traits (SvdTraits()).
     geom_traits(kernel));
}

// deprecated API
template <typename ConcurrencyTag,
	  typename InputIterator,
          typename PointMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::jet_smooth_point_set(), please update your code")
void
jet_smooth_point_set(
  InputIterator first, ///< iterator over the first input point
  InputIterator beyond, ///< past-the-end iterator
  PointMap point_map, ///< property map: value_type of InputIterator -> Point_3
  unsigned int k, ///< number of neighbors.
  const unsigned int degree_fitting = 2,
  const unsigned int degree_monge = 2)
{
  CGAL::Iterator_range<InputIterator> points (first, beyond);
  return jet_smooth_point_set<ConcurrencyTag>
    (points,
     k,
     CGAL::parameters::point_map (point_map).
     degree_fitting (degree_fitting).
     degree_monge (degree_monge));
}

// deprecated API
template <typename ConcurrencyTag,
	  typename InputIterator
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::jet_smooth_point_set(), please update your code")
void
jet_smooth_point_set(
  InputIterator first, ///< iterator over the first input point
  InputIterator beyond, ///< past-the-end iterator
  unsigned int k, ///< number of neighbors.
  const unsigned int degree_fitting = 2,
  const unsigned int degree_monge = 2)
{
  CGAL::Iterator_range<InputIterator> points (first, beyond);
  return jet_smooth_point_set<ConcurrencyTag>
    (points,
     k,
     CGAL::parameters::degree_fitting (degree_fitting).
     degree_monge (degree_monge));
}
#endif // CGAL Eigen / Lapack
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond
  

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_JET_SMOOTH_POINT_SET_H
