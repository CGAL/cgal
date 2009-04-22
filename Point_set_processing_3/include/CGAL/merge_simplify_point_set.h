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
// Author(s) : Laurent Saboret and Nader Salman

#ifndef CGAL_MERGE_SIMPLIFY_POINT_SET_H
#define CGAL_MERGE_SIMPLIFY_POINT_SET_H

#include <iterator>
#include <set>
#include <deque>
#include <algorithm>
#include <cmath>

CGAL_BEGIN_NAMESPACE


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace CGALi {


/// Utility class for merge_simplify_point_set():
/// Less_epsilon_points_3 defines a 3D points order / 2 points are equal
/// iff they belong to the same cell of a grid of cell size = epsilon.
template <class Point_3>
struct Less_epsilon_points_3
{
private:

    double m_epsilon;

public:

    Less_epsilon_points_3 (double epsilon)
        : m_epsilon (epsilon)
    {
        CGAL_precondition(epsilon > 0);
    }

    bool operator() (const Point_3& a, const Point_3& b) const
    {
        // Round points to multiples of m_epsilon, then compare.
        Point_3 rounded_a(round_epsilon(a.x(), m_epsilon),
                          round_epsilon(a.y(), m_epsilon),
                          round_epsilon(a.z(), m_epsilon));
        Point_3 rounded_b(round_epsilon(b.x(), m_epsilon),
                          round_epsilon(b.y(), m_epsilon),
                          round_epsilon(b.z(), m_epsilon));
        return (rounded_a < rounded_b);
    }

private:

    // Round number to multiples of epsilon
    static inline double round_epsilon(double value, double epsilon)
    {
        return std::floor(value/epsilon) * epsilon;
    }
};


} /* namespace CGALi */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Utility class for merge_simplify_point_set():
/// 3D points set which allows at most 1 point per cell
/// of a grid of cell size = epsilon.
///
/// Warning:
/// This class is a container sorted wrt points position
/// => you should not modify directly the order or the position of points.

template <class Point_3>
class Epsilon_point_set_3
  : public std::set<Point_3, CGALi::Less_epsilon_points_3<Point_3> >
{
private:

    // superclass
    typedef std::set<Point_3, CGALi::Less_epsilon_points_3<Point_3> > Base;

public:

    Epsilon_point_set_3 (double epsilon)
        : Base( CGALi::Less_epsilon_points_3<Point_3>(epsilon) )
    {
        CGAL_precondition(epsilon > 0);
    }

    // default copy constructor, operator =() and destructor are fine.
};


/// Merge points which belong to the same cell of a grid of cell size = epsilon.
/// This variant requires the kernel.
///
/// @commentheading Precondition: epsilon > 0.
///
/// @commentheading Template Parameters:
/// @param InputIterator value_type must be convertible to Point_3<Kernel>.
/// @param OutputIterator value_type must be convertible from InputIterator's value_type.
/// @param Kernel Geometric traits class.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator,
          typename Kernel
>
OutputIterator
merge_simplify_point_set(
          InputIterator first,      ///< iterator over the first input point.
          InputIterator beyond,     ///< past-the-end iterator over input points.
          OutputIterator output,    ///< iterator over the first output point.
          double epsilon,           ///< tolerance value when comparing 3D points.
          const Kernel& kernel)     ///< geometric traits.
{
    typedef typename std::iterator_traits<InputIterator>::value_type Point;

    CGAL_precondition(epsilon > 0);

    // Merge points which belong to the same cell of a grid of cell size = epsilon.
    Epsilon_point_set_3<Point> points_to_keep(epsilon);
    points_to_keep.insert(first, beyond);

    // Copy merged points to output
    output = std::copy(points_to_keep.begin(), points_to_keep.end(), output);
    return output;
}

/// Merge points which belong to the same cell of a grid of cell size = epsilon.
/// This variant is mutating the input point set and requires the kernel.
///
/// Warning:
/// This method modifies the order of points, thus
/// should not be called on sorted containers.
///
/// @commentheading Precondition: epsilon > 0.
///
/// @commentheading Template Parameters:
/// @param ForwardIterator value_type must be convertible to Point_3<Kernel>.
/// @param Kernel Geometric traits class.
///
/// @return First iterator to remove (see erase-remove idiom).
template <typename ForwardIterator,
          typename Kernel
>
ForwardIterator
merge_simplify_point_set(
           ForwardIterator first,   ///< iterator over the first input/output point.
           ForwardIterator beyond,  ///< past-the-end iterator.
           double epsilon,          ///< tolerance value when comparing 3D points.
           const Kernel& kernel)    ///< geometric traits.

{
    typedef typename std::iterator_traits<ForwardIterator>::value_type Point;

    CGAL_precondition(epsilon > 0);

    // Merge points which belong to the same cell of a grid of cell size = epsilon.
    // points_to_keep will contain 1 point per cell; the others will be in points_to_remove.
    Epsilon_point_set_3<Point> points_to_keep(epsilon);
    std::deque<Point> points_to_remove;
    for (ForwardIterator it=first ; it != beyond ; it++)
    {
        std::pair<Epsilon_point_set_3<Point>::iterator,bool> result;
        result = points_to_keep.insert(*it);
        if (!result.second) // if not inserted
          points_to_remove.push_back(*it);
    }

    // Replace [first, beyond) range by the content of points_to_keep, then points_to_remove.
    ForwardIterator first_iterator_to_remove =
    std::copy(points_to_keep.begin(), points_to_keep.end(), first);
    std::copy(points_to_remove.begin(), points_to_remove.end(), first_iterator_to_remove);

    return first_iterator_to_remove;
}

/// Merge points which belong to the same cell of a grid of cell size = epsilon.
/// This variant deduces the kernel from iterator types.
///
/// @commentheading Precondition: epsilon > 0.
///
/// @commentheading Template Parameters:
/// @param InputIterator value_type must be convertible to Point_3<Kernel>.
/// @param OutputIterator value_type must be convertible from InputIterator's value_type.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator
>
OutputIterator
merge_simplify_point_set(
           InputIterator first,   ///< iterator over the first input point.
           InputIterator beyond,  ///< past-the-end iterator over input points.
           OutputIterator output, ///< iterator over the first output point.
           double epsilon)        ///< tolerance value when comparing 3D points.
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return merge_simplify_point_set(first,beyond,output,epsilon,Kernel());
}

/// Merge points which belong to the same cell of a grid of cell size = epsilon.
/// This variant is mutating the input point set and deduces the kernel from iterator types.
///
/// Warning:
/// This method modifies the order of points, thus
/// should not be called on sorted containers.
///
/// @commentheading Precondition: epsilon > 0.
///
/// @commentheading Template Parameters:
/// @param ForwardIterator value_type must be convertible to Point_3<Kernel>.
///
/// @return First iterator to remove (see erase-remove idiom).
template <typename ForwardIterator>
ForwardIterator
merge_simplify_point_set(
       ForwardIterator first,       ///< iterator over the first input/output point
       ForwardIterator beyond,      ///< past-the-end iterator
       double epsilon)              ///< tolerance value when comparing 3D points
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return merge_simplify_point_set(first,beyond,epsilon,Kernel());
}


CGAL_END_NAMESPACE

#endif // CGAL_MERGE_SIMPLIFY_POINT_SET_H

