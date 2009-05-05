// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute point_it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s) : Nader Salman

#ifndef CGAL_REMOVE_OUTLIERS_WRT_MEDIAN_KNN_SQ_DISTANCE_H
#define CGAL_REMOVE_OUTLIERS_WRT_MEDIAN_KNN_SQ_DISTANCE_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <iterator>
#include <algorithm>
#include <map>

CGAL_BEGIN_NAMESPACE


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace CGALi {


/// Utility function for remove_outliers_wrt_median_knn_sq_distance():
/// Compute median squared distance to the K nearest neighbors.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param Kernel Geometric traits class.
/// @param Tree KD-tree.
///
/// @return computed distance.
template < typename Kernel,
           typename Tree >
typename Kernel::FT
compute_median_knn_sq_distance_3(
                     const typename Kernel::Point_3& query, ///< 3D point to project
                     Tree& tree,                            ///< KD-tree
                     unsigned int k)                        ///< number of neighbors
{
    // geometric types
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;

    // types for K nearest neighbors search
    typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;

    // Gather set of (k+1) neighboring points.
    // Perform k+1 queries (if in point set, the query point is
    // output first). Search may be aborted if k is greater
    // than number of input points.
    std::vector<Point> points; points.reserve(k+1);
    Neighbor_search search(tree,query,k+1);
    Search_iterator search_iterator = search.begin();
    unsigned int i;
    unsigned int k_median = k/2;
    for(i=0;i<(k_median+1);i++)/*for(i=0;i<(k+1);i++)*/
    {
        if(search_iterator == search.end())
            break; // premature ending
        points.push_back(search_iterator->first);
        search_iterator++;
    }
    CGAL_precondition(points.size() >= 1);

    // compute median squared distance
    typename Kernel::Compute_squared_distance_3 sqd;
    FT sq_distance = 0;
    //for(typename std::vector<Point>::iterator neighbor = points.begin(); neighbor != points.end(); neighbor++)
       // sq_distance += sqd(*neighbor, query);
    //sq_distance /= FT(points.size());

    //if((points.size() & 1) == 1)
        return sqd(points[k_median]/*[points.size()/2]*/, query);
    //else
    //    return (sqd(points[points.size()/2],query)+ sqd(points[(points.size()/2)+1], query))/2;
}


} /* namespace CGALi */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Remove outliers:
/// - compute median squared distance to the K nearest neighbors,
/// - sort the points in increasing order of computed distance.
///
/// This method modifies the order of input points, and returns 
/// an iterator over the first point to remove (see erase-remove idiom).
/// Warning: this method should not be called on sorted containers.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param ForwardIterator value_type must be convertible to Point_3<Kernel>.
/// @param Kernel Geometric traits class. It can be omitted and deduced automatically from the iterator type.
///
/// @return iterator over the first point to remove.

// This variant requires the kernel.
template <typename ForwardIterator,
          typename Kernel
>
ForwardIterator
remove_outliers_wrt_median_knn_sq_distance(
                     ForwardIterator first,     ///< iterator over the first input/output point.
                     ForwardIterator beyond,    ///< past-the-end iterator.
                     unsigned int k,            ///< number of neighbors.
                     const Kernel& kernel,      ///< geometric traits.
                     double threshold_percent)  ///< percentage of points to remove.
{
    // geometric types
    typedef typename Kernel::FT FT;
    typedef typename std::iterator_traits<ForwardIterator>::value_type Point;

    // types for K nearest neighbors search structure
    typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::Tree Tree;
    typedef typename Neighbor_search::iterator Search_iterator;

    // precondition: at least one element in the container.
    // to fix: should have at least three distinct points
    // but this is costly to check
    CGAL_precondition(first != beyond);

    // precondition: at least 2 nearest neighbors
    CGAL_precondition(k >= 2);

    CGAL_precondition(threshold_percent >= 0 && threshold_percent <= 100);

    // instanciate a KD-tree search
    Tree tree(first,beyond);

    // iterate over input points and add them to multimap sorted by distance to k
    std::multimap<FT,Point> sorted_points;
    for(ForwardIterator point_it = first; point_it != beyond; point_it++)
    {
      FT sq_distance = CGALi::compute_median_knn_sq_distance_3<Kernel>(*point_it,tree,k);
      sorted_points.insert( std::make_pair(sq_distance,*point_it) );
    }

    // Replace [first, beyond) range by the multimap content.
    // Return the iterator after the (100-threshold_percent) % best points.
    ForwardIterator first_point_to_remove = beyond;
    ForwardIterator dst = first;
    int first_index_to_remove = int(double(sorted_points.size()) * ((100.0-threshold_percent)/100.0));
    typename std::multimap<FT,Point>::iterator src;
    int index;
    for (src = sorted_points.begin(), index = 0;
         src != sorted_points.end();
         ++src, ++index)
    {
      *dst++ = src->second;
      if (index == first_index_to_remove)
        first_point_to_remove = dst;
    }

    return first_point_to_remove;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator>
ForwardIterator
remove_outliers_wrt_median_knn_sq_distance(
                     ForwardIterator first,     ///< iterator over the first input/output point.
                     ForwardIterator beyond,    ///< past-the-end iterator.
                     unsigned int k,            ///< number of neighbors.
                     double threshold_percent)  ///< percentage of points to remove.
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return remove_outliers_wrt_median_knn_sq_distance(first,beyond,k,Kernel(),threshold_percent);
}
/// @endcond


CGAL_END_NAMESPACE

#endif // CGAL_REMOVE_OUTLIERS_WRT_MEDIAN_KNN_SQ_DISTANCE_H

