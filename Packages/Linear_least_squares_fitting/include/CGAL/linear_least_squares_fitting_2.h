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
// Author(s)     : Pierre Alliez and Sylvain Pion

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_2
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_2

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <iterator>
#include <list>
#include <string>
#include <CGAL/eigen_2.h>
#include <CGAL/centroid.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// fits a line to a 2D point set
// returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best (zero variance orthogonally to the fitting line)
//  0 is worst (isotropic case, returns a line with default direction)
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const K& k,                 // kernel
                               const typename K::Point_2*)
{
  typedef typename K::FT       FT;
  typedef typename K::Point_2  Point;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Line_2   Line;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K());

  // assemble covariance matrix as a
  // semi-definite matrix. 
  // Matrix numbering:
  // 0
  // 1 2
  FT covariance[3] = {0,0,0};
  std::pair<FT,FT> eigen_values;
  std::pair<Vector,Vector> eigen_vectors;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Point& p = *it;
    Vector d = p - c;
    covariance[0] += d.x() * d.x();
    covariance[1] += d.x() * d.y();
    covariance[2] += d.y() * d.y();
  }

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  eigen_symmetric_2<K>(covariance,eigen_vectors,eigen_values);

  // assert eigen values are positives or null
  CGAL_assertion(eigen_values.first >= 0 && 
                 eigen_values.second >= 0);

  // check unicity and build fitting line accordingly
  if(eigen_values.first != eigen_values.second)
  {
    // regular case
    line = Line(c,eigen_vectors.first);
    return (FT)1.0 - eigen_values.second / eigen_values.first;
  } // end regular case
  else
  {
    // isotropic case (infinite number of directions)
    // by default: assemble a line that goes through 
    // the centroid and with a default horizontal vector.
    line = Line(c,Vector(1,0));
    return (FT)0.0;
  } // end isotropic case
  // end linear_least_squares_fitting_2 for point set

} // namespace CGALi


template < typename InputIterator, 
           typename K >
inline
typename K::FT
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               typename K::Line_2& line,
                               typename K::Point_2& centroid,
                               const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return CGALi::linear_least_squares_fitting_2(begin, end, line,
                                               centroid, k, (Value_type*) NULL);
}


template < typename InputIterator, 
           typename K >
inline
typename K::FT
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               typename K::Line_2& line,
                               const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typename K::Point_2 centroid;
  return CGALi::linear_least_squares_fitting_2(begin, end, line,
                                               centroid, k,(Value_type*) NULL);
}


// these ones deduce the kernel from the 
// points in container

template < typename InputIterator, 
           typename Line,
           typename Point>
inline
typename Kernel_traits<Line>::Kernel::FT
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               Line& line,
                               Point& centroid)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel K;
  return CGAL::linear_least_squares_fitting_2(begin,end,line,centroid,K());
}

template < typename InputIterator, 
           typename Line >
inline
typename Kernel_traits<Line>::Kernel::FT
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               Line& line)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel K;
  return CGAL::linear_least_squares_fitting_2(begin,end,line,K());
}

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_2


