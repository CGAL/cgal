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

#ifndef CGAL_MERGE_SIMPLIFY_NEAREST_POINT_SET_H
#define CGAL_MERGE_SIMPLIFY_NEAREST_POINT_SET_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <algorithm>
#include <map>
#include <set>
#include <deque>
#include <cmath>

CGAL_BEGIN_NAMESPACE

template <typename InputIterator,
          typename OutputIterator,
          typename Kernel
>
OutputIterator
merge_simplify_nearest_point_set(
          InputIterator first,      ///< iterator over the first input point.
          InputIterator beyond,     ///< past-the-end iterator.
          OutputIterator output,    ///< iterator over the first output point.
          unsigned int k,           ///< number of neighbors.
          double epsilon,           ///< tolerance value when comparing 3D points.
          const Kernel& kernel)     ///< geometric traits.
{
    // geometric types
    typedef typename Kernel::FT FT;
    typedef typename std::iterator_traits<InputIterator>::value_type Point;


    // types for K nearest neighbors search
    typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
    typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
    typedef typename Neighbor_search::Tree Tree;
    typedef typename Neighbor_search::iterator Search_iterator;

    // preconditions
    CGAL_point_set_processing_precondition(first != beyond);
    CGAL_point_set_processing_precondition(k >= 1);
    CGAL_point_set_processing_precondition(epsilon > 0);

    // instanciate a KD-tree search
    Tree tree(first,beyond);

    // Merge points which belong to the same cell of a grid of cell size = epsilon.
    // points_to_keep will contain 1 point per cell; the others will be in points_to_remove.
    std::set<Point>   points_to_keep;

    typename Kernel::Compute_squared_distance_3 sqd;

    for (InputIterator it=first ; it != beyond ; it++)
    {
        std::vector<Point> points; points.reserve(k+1);
        Neighbor_search search(tree,*it,k+1);
        Search_iterator search_iterator = search.begin();
        unsigned int i;

        for(i=0;i<(k+1);i++)
        {
            if(search_iterator == search.end())
                break; // premature ending
            if(sqd(search_iterator->first, *it)<= epsilon)
            {
                (*it).merge(search_iterator->first);
                search_iterator++;
            }

            else
            {
                points_to_keep.insert(search_iterator->first);
                search_iterator++;
            }
        }
        points_to_keep.insert(*it);
    }

    // Copy merged points to output
    output = std::copy(points_to_keep.begin(), points_to_keep.end(), output);
    return output;
}



template <typename InputIterator,
          typename OutputIterator
>
OutputIterator
merge_simplify_nearest_point_set(
          InputIterator first,      ///< iterator over the first input point
          InputIterator beyond,     ///< past-the-end iterator
          OutputIterator output,    ///< iterator over the first output point
          unsigned int k,           ///< number of neighbors
          double epsilon)           ///< tolerance value when comparing 3D points
{
    typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
    typedef typename CGAL::Kernel_traits<Value_type>::Kernel Kernel;
    return merge_simplify_nearest_point_set(first,beyond,output,k,epsilon,Kernel());
}


CGAL_END_NAMESPACE

#endif // CGAL_MERGE_SIMPLIFY_NEAREST_POINT_SET_H

