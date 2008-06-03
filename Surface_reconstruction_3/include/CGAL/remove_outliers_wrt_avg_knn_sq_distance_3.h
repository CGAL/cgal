// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
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
// Author(s) : Laurent Saboret and Nader Salman

#ifndef CGAL_REMOVE_OUTLIERS_WRT_AVG_KNN_SQ_DISTANCE_3_H
#define CGAL_REMOVE_OUTLIERS_WRT_AVG_KNN_SQ_DISTANCE_3_H

#include <CGAL/basic.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <iterator>
#include <map>

CGAL_BEGIN_NAMESPACE


/// Compute average squared distance to the K nearest neighbors.
///
/// Precondition: KNN >= 2.
///
/// @heading Parameters:
/// @param Kernel Geometric traits class.
/// @param Tree KD-tree.
///
/// @return computed distance.
template < typename Kernel,
           typename Tree >
typename Kernel::FT
compute_avg_knn_sq_distance_3(
                     const typename Kernel::Point_3& query, ///< 3D point to project
                     Tree& tree,                            ///< KD-tree
                     unsigned int KNN)                      ///< number of neighbors
{
    // geometric types
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;

    // types for K nearest neighbors search
    typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::iterator Search_iterator;

    // gather set of (KNN+1) neighboring points
    // performs KNN + 1 queries (if unique the query point is
    // output first). search may be aborted when KNN is greater
    // than number of input points
    std::vector<Point> points;
    Neighbor_search search(tree,query,KNN+1);
    Search_iterator search_iterator = search.begin();
    unsigned int i;
    for(i=0;i<(KNN+1);i++)
    {
        if(search_iterator == search.end())
            break; // premature ending
        points.push_back(search_iterator->first);
        search_iterator++;
    }
    CGAL_precondition(points.size() >= 1);

    // compute squared distance
    typename Kernel::Compute_squared_distance_3 sqd;
    FT sq_distance = 0;
    for(typename std::vector<Point>::iterator neighbor = points.begin(); neighbor != points.end(); neighbor++)
        sq_distance += sqd(*neighbor, query);
    sq_distance /= FT(KNN);
    return sq_distance;
}

/// Remove outliers:
/// - compute average squared distance to the K nearest neighbors,
/// - percentage of points to remove.
/// This variant requires the kernel.
///
/// Precondition: KNN >= 2.
///
/// @heading Parameters:
/// @param InputIterator value_type is Point_3.
/// @param OutputIterator value_type is Point_3.
/// @param Kernel Geometric traits class.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator,
          typename Kernel
>
OutputIterator
remove_outliers_wrt_avg_knn_sq_distance_3(
                     InputIterator first,       ///< input points
                     InputIterator beyond,
                     OutputIterator output,     ///< output points
                     unsigned int KNN,          ///< number of neighbors
                     const Kernel& /*kernel*/,
                     double threshold_percent)  ///< percentage of points to remove
{
  // geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;

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
  CGAL_precondition(KNN >= 2);

  CGAL_precondition(threshold_percent >= 0 && threshold_percent <= 100);

  // instanciate a KD-tree search
  Tree tree(first,beyond);

  // iterate over input points and add them to multimap sorted by distance to knn
  std::multimap<FT,InputIterator> map;
  for(InputIterator point_it = first; point_it != beyond; point_it++)
  {
    FT sq_distance = compute_avg_knn_sq_distance_3<Kernel,Tree>(*point_it,tree,KNN);
    map.insert( std::pair<FT,InputIterator>(sq_distance,point_it) );
  }

  // output (100-threshold_percent) % best points
  typename std::multimap<FT,InputIterator>::iterator map_it;
  int index;
  int last = int(double(map.size()) * ((100.0-threshold_percent)/100.0));
  for(map_it = map.begin(), index=0; index < last; ++map_it, ++index)
  {
    InputIterator point_it = map_it->second;
    *output = *point_it;
    output++;
  }
  return output;
}

/// Remove outliers:
/// - compute average squared distance to the K nearest neighbors,
/// - percentage of points to remove.
/// This function is mutating the input point set.
/// This variant requires the kernel.
///
/// Precondition: KNN >= 2.
///
/// @heading Parameters:
/// @param ForwardIterator value_type is Point_3.
/// @param Kernel Geometric traits class.
template <typename ForwardIterator,
          typename Kernel
>
void
remove_outliers_wrt_avg_knn_sq_distance_3(
                     ForwardIterator first,     ///< input/output points
                     ForwardIterator beyond,
                     unsigned int KNN,          ///< number of neighbors
                     const Kernel& /*kernel*/,
                     double threshold_percent)  ///< percentage of points to remove
{
  CGAL_precondition(false); // nyi
}


/// Remove outliers:
/// - compute average squared distance to the K nearest neighbors,
/// - percentage of points to remove.
/// This variant deduces the kernel from iterator types.
///
/// Precondition: KNN >= 2.
///
/// @heading Parameters:
/// @param InputIterator value_type is Point_3.
/// @param OutputIterator value_type is Point_3.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator
>
OutputIterator
remove_outliers_wrt_avg_knn_sq_distance_3(
                     InputIterator first,       ///< input points
                     InputIterator beyond,
                     OutputIterator output,     ///< output points
                     unsigned int KNN,          ///< number of neighbors
                     double threshold_percent)  ///< percentage of points to remove
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return remove_outliers_wrt_avg_knn_sq_distance_3(first,beyond,output,KNN,Kernel(),threshold_percent);
}

/// Remove outliers:
/// - compute average squared distance to the K nearest neighbors,
/// - percentage of points to remove.
/// This function is mutating the input point set.
/// This variant deduces the kernel from iterator types.
///
/// Precondition: KNN >= 2.
///
/// @heading Parameters:
/// @param ForwardIterator value_type is Point_3.
template <typename ForwardIterator>
void
remove_outliers_wrt_avg_knn_sq_distance_3(
                     ForwardIterator first,     ///< input/output points
                     ForwardIterator beyond,
                     unsigned int KNN,          ///< number of neighbors
                     double threshold_percent)  ///< percentage of points to remove
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  remove_outliers_wrt_avg_knn_sq_distance_3(first,beyond,KNN,Kernel(),threshold_percent);
}


CGAL_END_NAMESPACE

#endif // CGAL_REMOVE_OUTLIERS_WRT_AVG_KNN_SQ_DISTANCE_3_H

