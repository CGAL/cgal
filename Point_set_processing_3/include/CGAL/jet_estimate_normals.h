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
// Author(s) : Pierre Alliez and Laurent Saboret and Marc Pouget and Frederic Cazals

#ifndef CGAL_JET_ESTIMATE_NORMALS_H
#define CGAL_JET_ESTIMATE_NORMALS_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/trace.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Memory_sizer.h>

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


/// Estimates normal direction using jet fitting
/// on the k nearest neighbors.
///
/// \pre `k >= 2`
///
/// @tparam Kernel Geometric traits class.
/// @tparam Tree KD-tree.
///
/// @return Computed normal. Orientation is random.
template < typename Kernel,
           typename SvdTraits,
           typename Tree
>
typename Kernel::Vector_3
jet_estimate_normal(const typename Kernel::Point_3& query, ///< point to compute the normal at
                    Tree& tree, ///< KD-tree
                    unsigned int k, ///< number of neighbors
                    unsigned int degree_fitting)
{
  // basic geometric types
  typedef typename Kernel::Point_3  Point;

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
  // Perform k+1 queries (as in point set, the query point is
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
  const unsigned int degree_monge = 1; // we seek for normal and not more.
  Monge_form monge_form = monge_fit(points.begin(), points.end(),
                                    degree_fitting, degree_monge);

  // output normal vector (already normalized in monge form)
  return monge_form.normal_direction();
}

#ifdef CGAL_LINKED_WITH_TBB
  template <typename Kernel, typename SvdTraits, typename Tree>
  class Jet_estimate_normals {
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    const Tree& tree;
    const unsigned int k;
    const unsigned int degree_fitting;
    const std::vector<Point>& input;
    std::vector<Vector>& output;

  public:
    Jet_estimate_normals(Tree& tree, unsigned int k, std::vector<Point>& points,
		     unsigned int degree_fitting, std::vector<Vector>& output)
      : tree(tree), k (k), degree_fitting (degree_fitting), input (points), output (output)
    { }
    
    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for( std::size_t i = r.begin(); i != r.end(); ++i)
	output[i] = CGAL::internal::jet_estimate_normal<Kernel,SvdTraits>(input[i], tree, k, degree_fitting);
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
   Estimates normal directions of the range of `points`
   using jet fitting on the k nearest neighbors.
   The output normals are randomly oriented.

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
     \cgalParamBegin{normal_map} a model of `ReadWritePropertyMap` with value type
     `geom_traits::Vector_3`.\cgalParamEnd
     \cgalParamBegin{degree_fitting} degree of jet fitting.\cgalParamEnd
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
jet_estimate_normals(
  PointRange& points,
  unsigned int k,
  const NamedParameters& np)
{
  using boost::choose_param;
  
  CGAL_TRACE("Calls jet_estimate_normals()\n");

  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::type NormalMap;
  typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
  typedef typename GetSvdTraits<NamedParameters>::type SvdTraits;

  CGAL_static_assertion_msg(!(boost::is_same<NormalMap,
                              typename Point_set_processing_3::GetNormalMap<PointRange, NamedParameters>::NoMap>::value),
                            "Error: no normal map");
  CGAL_static_assertion_msg(!(boost::is_same<SvdTraits,
                              typename GetSvdTraits<NamedParameters>::NoTraits>::value),
                            "Error: no SVD traits");

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  NormalMap normal_map = choose_param(get_param(np, internal_np::normal_map), NormalMap());
  unsigned int degree_fitting = choose_param(get_param(np, internal_np::degree_fitting), 2);

  typedef typename Kernel::Point_3 Point;

  // Input points types
  typedef typename boost::property_traits<NormalMap>::value_type Vector;

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

  std::size_t memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("  Creates KD-tree\n");

  typename PointRange::iterator it;

  // Instanciate a KD-tree search.
  // Note: We have to convert each input iterator to Point_3.
  std::vector<Point> kd_tree_points; 
  for(it = points.begin(); it != points.end(); it++)
    kd_tree_points.push_back(get(point_map, *it));
  Tree tree(kd_tree_points.begin(), kd_tree_points.end());

  memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("  Computes normals\n");

  // iterate over input points, compute and output normal
  // vectors (already normalized)
#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
			     "Parallel_tag is enabled but TBB is unavailable.");
#else
   if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
   {
     std::vector<Vector> normals (kd_tree_points.size ());
     CGAL::internal::Jet_estimate_normals<Kernel, SvdTraits, Tree>
       f (tree, k, kd_tree_points, degree_fitting, normals);
     tbb::parallel_for(tbb::blocked_range<size_t>(0, kd_tree_points.size ()), f);
     unsigned int i = 0;
     for(it = points.begin(); it != points.end(); ++ it, ++ i)
       {
	 put (normal_map, *it, normals[i]);
       }
   }
   else
#endif
     {
       for(it = points.begin(); it != points.end(); it++)
	 {
	   Vector normal = internal::jet_estimate_normal<Kernel,SvdTraits,Tree>(
										get(point_map,*it), 
										tree, k, degree_fitting);

	   put(normal_map, *it, normal); // normal_map[it] = normal
    
	 }
     }


  memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("End of jet_estimate_normals()\n");
}


/// \cond SKIP_IN_MANUAL
// variant with default NP
template <typename ConcurrencyTag,
	  typename PointRange>
void
jet_estimate_normals(
  PointRange& points,
  unsigned int k) ///< number of neighbors.
{
  jet_estimate_normals<ConcurrencyTag>
    (points, k, CGAL::Point_set_processing_3::parameters::all_default(points));
}

#ifndef CGAL_NO_DEPRECATED_CODE
// deprecated API
template <typename ConcurrencyTag,
	  typename ForwardIterator,
          typename PointMap,
          typename NormalMap,
          typename Kernel,
          typename SvdTraits
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::jet_estimate_normals(), please update your code")
void
jet_estimate_normals(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
  NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3.
  unsigned int k, ///< number of neighbors.
  const Kernel& /*kernel*/, ///< geometric traits.
  unsigned int degree_fitting = 2) ///< fitting degree
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return jet_estimate_normals<ConcurrencyTag>
    (points,
     k,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map).
     degree_fitting (degree_fitting).
     geom_traits(Kernel()));
}
  
#if defined(CGAL_EIGEN3_ENABLED) || defined(CGAL_LAPACK_ENABLED)
// deprecated API
template <typename ConcurrencyTag,
	  typename ForwardIterator,
          typename PointMap,
          typename NormalMap,
          typename Kernel
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::jet_estimate_normals(), please update your code")
void
jet_estimate_normals(
  ForwardIterator first,
  ForwardIterator beyond,
  PointMap point_map,
  NormalMap normal_map,
  unsigned int k,
  const Kernel& kernel,
  unsigned int degree_fitting = 2)
{
  #ifdef CGAL_EIGEN3_ENABLED
  typedef Eigen_svd SvdTraits;
  #else
  typedef Lapack_svd SvdTraits;
  #endif

  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return jet_estimate_normals<ConcurrencyTag>
    (points,
     k,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map).
     degree_fitting (degree_fitting).
     svd_traits (SvdTraits()).
     geom_traits(kernel));
}

// deprecated API
template <typename ConcurrencyTag,
	  typename ForwardIterator,
          typename PointMap,
          typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::jet_estimate_normals(), please update your code")
void
jet_estimate_normals(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointMap point_map, ///< property map: value_type of ForwardIterator -> Point_3.
  NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3.
  unsigned int k, ///< number of neighbors.
  unsigned int degree_fitting = 2)
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return jet_estimate_normals<ConcurrencyTag>
    (points,
     k,
     CGAL::parameters::point_map (point_map).
     normal_map (normal_map).
     degree_fitting (degree_fitting));
}
  
// deprecated API
template <typename ConcurrencyTag,
	  typename ForwardIterator,
          typename NormalMap
>
CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::jet_estimate_normals(), please update your code")
void
jet_estimate_normals(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  NormalMap normal_map, ///< property map: value_type of ForwardIterator -> Vector_3.
  unsigned int k, ///< number of neighbors.
  unsigned int degree_fitting = 2)
{
  CGAL::Iterator_range<ForwardIterator> points (first, beyond);
  return jet_estimate_normals<ConcurrencyTag>
    (points,
     k,
     CGAL::parameters::normal_map (normal_map).
     degree_fitting (degree_fitting));
}
#endif // CGAL Eigen / Lapack
#endif // CGAL_NO_DEPRECATED_CODE
/// \endcond

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_JET_ESTIMATE_NORMALS_H
