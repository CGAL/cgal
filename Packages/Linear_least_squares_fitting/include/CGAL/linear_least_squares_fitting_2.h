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
// $Source:
// /CVSROOT/CGAL/Packages/Interpolation/include/CGAL/interpolation_functions.h,v $
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Pierre Alliez

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_2
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_2

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <iterator>
#include <list>
#include <string>
#include <CGAL/eigen.h>
#include <CGAL/centroid.h>

CGAL_BEGIN_NAMESPACE

// fit a line to a 2D point set
// return a fitting quality (1 - lambda_min/lambda_max):
//  1 is best (zero variance orthogonally to the fitting line)
//  0 is worst (isotropic case, return a line with default orientation)
template < typename InputIterator, 
           typename Line, 
           typename K >
typename K::FT
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               Line& line,
                               const K& k)
{
  typedef typename K::FT          FT;
  typedef typename K::Point_2     Point;
  typedef typename K::Vector_2    Vector;
  typedef typename K::Segment_2   Segment;
  typedef typename K::Direction_2 Direction;

  // precondition: at least one element in the container.
  CGAL_precondition(begin != end);

  // compute centroid
  Point c = centroid(begin,end,K());

  // assemble covariance matrix as a
  // semi-definite matrix. 
  // Matrix numbering:
  // 0          
  // 1 2
  FT covariance[3] = {0,0,0};
  FT eigen_values[2] = {0,0};
  // these should be typed Vector_2
  // but it requires working more on eigen.h
  FT eigen_vectors[4] = {0,0,0,0};
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
  // TODO: use explicit formula instead.
  eigen_semi_definite_symmetric(covariance,2,eigen_vectors,eigen_values);

  // assert eigen values are positives
  CGAL_assertion(eigen_values[0] >= 0 && 
                       eigen_values[1] >= 0);

  // check unicity and build fitting line accordingly
  if(eigen_values[0] == eigen_values[1])
  {
    // regular case
    line = Line(c,Direction(eigen_vectors[0],eigen_vectors[1]));
    return (ft)1.0 - eigen_values[1] / eigen_values[0];
  } // end regular case
  else
  {
    // isotropic case (infinite number of directions)
    // by default: assemble a line that goes through 
    // the centroid and with a default direction.
    line = Line(c,Direction());
    return (FT)0.0;
  } // end isotropic case
} // end linear_least_squares_fitting_2


// this one deduces the kernel from the point in container
template < typename InputIterator, 
           typename Line >
typename Kernel_traits<Line>::Kernel::FT
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               Line& line)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel                    K;
  return CGAL::linear_least_squares_fitting_2(begin,end,line,K());
}

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_2


