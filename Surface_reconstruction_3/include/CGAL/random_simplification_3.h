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
//
// Author(s) : Laurent Saboret

#ifndef CGAL_RANDOM_SIMPLIFICATION_3_H
#define CGAL_RANDOM_SIMPLIFICATION_3_H

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>

CGAL_BEGIN_NAMESPACE


/// Delete random points.
/// This variant requires the kernel.
///
/// @heading Parameters:
/// @param InputIterator value_type must be convertible to OutputIterator's value_type.
/// @param OutputIterator value_type must be convertible to Point_3.
/// @param Kernel Geometric traits class.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator,
          typename Kernel
>
OutputIterator
random_simplification_3(
          InputIterator first,      ///< input points
          InputIterator beyond,
          OutputIterator output,    ///< output points
          double threshold_percent, ///< percentage of points to remove
          const Kernel& /*kernel*/)          
{
    typedef typename std::iterator_traits<InputIterator>::value_type Point;

    CGAL_precondition(threshold_percent >= 0 && threshold_percent <= 100);

    // Random shuffle
    std::vector<Point> points (first, beyond);
    std::random_shuffle (points.begin(), points.end());

    // Compute first iterator to remove
    int nb_points = points.size();
    int first_index_to_remove = int(double(nb_points) * ((100.0-threshold_percent)/100.0));

    // Output points up to first_index_to_remove
    output = std::copy(&points[0], &points[first_index_to_remove], output);
    return output;
}

/// Delete random points.
/// This function is mutating the input point set.
/// This variant requires the kernel.
///
/// Warning: 
/// This method modifies the order of points, thus
/// should not be called on sorted containers.
///
/// @heading Parameters:
/// @param ForwardIterator value_type must be convertible to Point_3.
/// @param Kernel Geometric traits class.
///
/// @return First iterator to remove (see erase-remove idiom).
template <typename ForwardIterator,
          typename Kernel
>
ForwardIterator
random_simplification_3(
           ForwardIterator first,     ///< input/output points
           ForwardIterator beyond,
           double threshold_percent,  ///< percentage of points to remove
           const Kernel& /*kernel*/)
        
{
    typedef typename std::iterator_traits<ForwardIterator>::value_type Point;

    CGAL_precondition(threshold_percent >= 0 && threshold_percent <= 100);

    // Random shuffle
    std::random_shuffle (first, beyond);

    // Compute first iterator to remove
    int nb_points = std::distance(first, beyond);
    int first_index_to_remove = int(double(nb_points) * ((100.0-threshold_percent)/100.0));
    ForwardIterator first_iterator_to_remove = first;
    std::advance(first_iterator_to_remove, first_index_to_remove);
    
    return first_iterator_to_remove;
}

/// Delete random points.
/// This variant deduces the kernel from iterator types.
///
/// @heading Parameters:
/// @param InputIterator value_type must be convertible to OutputIterator's value_type.
/// @param OutputIterator value_type must be convertible to Point_3.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator
>
OutputIterator
random_simplification_3(
           InputIterator first,       ///< input points
           InputIterator beyond,
           OutputIterator output,     ///< output points
           double threshold_percent)  ///< percentage of points to remove
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return random_simplification_3(first,beyond,output,threshold_percent,Kernel());
}

/// Delete random points.
/// This function is mutating the input point set.
/// This variant deduces the kernel from iterator types.
///
/// Warning: 
/// This method modifies the order of points, thus
/// should not be called on sorted containers.
///
/// @heading Parameters:
/// @param ForwardIterator value_type must be convertible to Point_3.
///
/// @return First iterator to remove (see erase-remove idiom).
template <typename ForwardIterator>
ForwardIterator
random_simplification_3(
       ForwardIterator first,     ///< input/output points
       ForwardIterator beyond,
       double threshold_percent)  ///< percentage of points to remove
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return random_simplification_3(first,beyond,threshold_percent,Kernel());
}


CGAL_END_NAMESPACE

#endif // CGAL_RANDOM_SIMPLIFICATION_3_H

