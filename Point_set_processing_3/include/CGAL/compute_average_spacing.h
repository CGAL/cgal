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
//
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_AVERAGE_SPACING_3_H
#define CGAL_AVERAGE_SPACING_3_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/Search_traits_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/assertions.h>

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

/// \ingroup PkgPointSetProcessingAlgorithms
/// Computes average spacing from k nearest neighbors.
///
/// \pre `k >= 2.`
///
/// @tparam Concurrency_tag enables sequential versus parallel algorithm.
///                         Possible values are `Sequential_tag`
///                         and `Parallel_tag`.
/// @tparam InputIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
///        It can be omitted if the value type of `InputIterator`  is convertible to `Point_3<Kernel>`.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointPMap`.
///
/// @return average spacing (scalar).

// This variant requires the kernel.
template <typename Concurrency_tag,
	  typename InputIterator,
          typename PointPMap,
          typename Kernel
>
typename Kernel::FT
compute_average_spacing(
  InputIterator first,  ///< iterator over the first input point.
  InputIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3
  unsigned int k, ///< number of neighbors.
  const Kernel& /*kernel*/) ///< geometric traits.
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point;

  // types for K nearest neighbors search structure
  typedef typename Kernel::FT FT;
  typedef Search_traits_3<Kernel> Tree_traits;
  typedef Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_point_set_processing_precondition(first != beyond);

  // precondition: at least 2 nearest neighbors
  CGAL_point_set_processing_precondition(k >= 2);

  // Instanciate a KD-tree search.
  // Note: We have to convert each input iterator to Point_3.
  std::vector<Point> kd_tree_points; 
  for(InputIterator it = first; it != beyond; it++)
    kd_tree_points.push_back(get(point_pmap, *it));
  Tree tree(kd_tree_points.begin(), kd_tree_points.end());

  // iterate over input points, compute and output normal
  // vectors (already normalized)
  FT sum_spacings = (FT)0.0;

#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<Concurrency_tag, Parallel_tag>::value),
			     "Parallel_tag is enabled but TBB is unavailable.");
#else
   if (boost::is_convertible<Concurrency_tag,Parallel_tag>::value)
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
       for(InputIterator it = first; it != beyond; it++)
	 {
	   sum_spacings += internal::compute_average_spacing<Kernel,Tree>(
									  get(point_pmap,*it),
									  tree,k);
	 }
     }
   
  // return average spacing
   return sum_spacings / (FT)(kd_tree_points.size ());
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename Concurrency_tag,
	  typename InputIterator,
          typename PointPMap
>
typename Kernel_traits<typename boost::property_traits<PointPMap>::value_type>::Kernel::FT
compute_average_spacing(
  InputIterator first,    ///< iterator over the first input point.
  InputIterator beyond,   ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3
  unsigned int k) ///< number of neighbors
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return compute_average_spacing<Concurrency_tag>(
    first,beyond,
    point_pmap,
    k,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
template < typename Concurrency_tag, typename InputIterator >
typename Kernel_traits<typename std::iterator_traits<InputIterator>::value_type>::Kernel::FT
compute_average_spacing(
  InputIterator first,    ///< iterator over the first input point.
  InputIterator beyond,   ///< past-the-end iterator over the input points.
  unsigned int k) ///< number of neighbors.
{
  return compute_average_spacing<Concurrency_tag>(
    first,beyond,
    make_identity_property_map(
    typename std::iterator_traits<InputIterator>::value_type()),
    k);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_AVERAGE_SPACING_3_H
