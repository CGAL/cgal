// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s) : Pierre Alliez and Sylvain Pion

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/centroid.h>
#include <CGAL/eigen_2.h>

#include <iterator>
#include <list>
#include <string>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// Fits a line to a 2D point set.
// Returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best  (zero variance orthogonally to the fitting line);
//  0 is worst (isotropic case, returns a line with horizontal
//              direction by default).
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const K&,                   // kernel
                               const typename K::Point_2*) // used for indirection
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Line_2   Line;
  typedef typename K::Point_2  Point;
  typedef typename K::Vector_2 Vector;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K());

  // assemble covariance matrix as a semi-definite matrix. 
  // Matrix numbering:
  // 0
  // 1 2
  FT covariance[3] = {0,0,0};
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Point& p = *it;
    Vector d = p - c; // centered data point
    covariance[0] += d.x() * d.x();
    covariance[1] += d.x() * d.y();
    covariance[2] += d.y() * d.y();
  }

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  std::pair<FT,FT> eigen_values;
  std::pair<Vector,Vector> eigen_vectors;
  CGALi::eigen_symmetric_2<K>(covariance, eigen_vectors, eigen_values);

  // check unicity and build fitting line accordingly
  if(eigen_values.first != eigen_values.second)
  {
    // regular case
    line = Line(c, eigen_vectors.first);
    return (FT)1.0 - eigen_values.second / eigen_values.first;
  } 
  else
  {
    // isotropic case (infinite number of directions)
    // by default: assemble a line that goes through 
    // the centroid and with a default horizontal vector.
    line = Line(c, Vector(1.0, 0.0));
    return (FT)0.0;
  } 
} // end linear_least_squares_fitting_2 for point set

} // end namespace CGALi


template < typename InputIterator, 
           typename K >
inline
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_2& line,
                               typename K::Point_2& centroid,
                               const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return CGALi::linear_least_squares_fitting_2(first, beyond, line,
                                               centroid, k, (Value_type*) NULL);
}

template < typename InputIterator, 
           typename K >
inline
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_2& line,
                               const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typename K::Point_2 centroid;
  return CGALi::linear_least_squares_fitting_2(first, beyond, line,
                                               centroid, k,(Value_type*) NULL);
}

// deduces the kernel from the points in container.
template < typename InputIterator, 
           typename Line,
           typename Point>
inline
typename Kernel_traits<Line>::Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               Line& line,
                               Point& centroid)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel K;
  return CGAL::linear_least_squares_fitting_2(first,beyond,line,centroid,K());
}

// does not return the centroid and deduces the kernel as well.
template < typename InputIterator, 
           typename Line >
inline
typename Kernel_traits<Line>::Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               Line& line)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel K;
  return CGAL::linear_least_squares_fitting_2(first,beyond,line,K());
}

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H
