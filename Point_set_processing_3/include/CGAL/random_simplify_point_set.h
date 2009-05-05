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

#ifndef CGAL_RANDOM_SIMPLIFY_POINT_SET_H
#define CGAL_RANDOM_SIMPLIFY_POINT_SET_H

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>

CGAL_BEGIN_NAMESPACE


/// Randomly deletes a user-specified fraction of the input points.
///
/// This method modifies the order of input points, and returns 
/// an iterator over the first point to remove (see erase-remove idiom).
/// Warning: this method should not be called on sorted containers.
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
random_simplify_point_set(
           ForwardIterator first,     ///< iterator over the first input/output point.
           ForwardIterator beyond,    ///< past-the-end iterator.
           double threshold_percent,  ///< percentage of points to remove.
           const Kernel& kernel)      ///< geometric traits.
{
    typedef typename std::iterator_traits<ForwardIterator>::value_type Point;

    CGAL_precondition(threshold_percent >= 0 && threshold_percent <= 100);

    // Random shuffle
    std::random_shuffle (first, beyond);

    // Compute first iterator to remove
    int nb_points = std::distance(first, beyond);
    int first_index_to_remove = int(double(nb_points) * ((100.0-threshold_percent)/100.0));
    ForwardIterator first_point_to_remove = first;
    std::advance(first_point_to_remove, first_index_to_remove);

    return first_point_to_remove;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator>
ForwardIterator
random_simplify_point_set(
       ForwardIterator first,     ///< iterator over the first input/output point
       ForwardIterator beyond,    ///< past-the-end iterator
       double threshold_percent)  ///< percentage of points to remove
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return random_simplify_point_set(first,beyond,threshold_percent,Kernel());
}
/// @endcond


CGAL_END_NAMESPACE

#endif // CGAL_RANDOM_SIMPLIFY_POINT_SET_H

