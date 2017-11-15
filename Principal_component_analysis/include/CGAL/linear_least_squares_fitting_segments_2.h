// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_2_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
#include <CGAL/centroid.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/PCA_util.h>

#include <iterator>
#include <list>
#include <cmath>

namespace CGAL {

namespace internal {
// Fits a line to a 2D segment set.
// Returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best  (zero variance orthogonally to the fitting line);
//  0 is worst (isotropic case, returns a line with horizontal
//              direction by default)

template < typename InputIterator, typename K, typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Segment_2*,// used for indirection
                               const K&,                   // kernel
			       const CGAL::Dimension_tag<1>& tag,
			       const DiagonalizeTraits&)
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Line_2   Line;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Segment_2 Segment;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K(),tag);
  // assemble covariance matrix as a semi-definite matrix. 
  // Matrix numbering:
  // 0
  // 1 2
  //Final combined covariance matrix for all segments and their combined mass
  FT mass = 0.0;
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0. }};

  // assemble 2nd order moment about the origin.  
  FT temp[4] = {1.0, 0.5, 0.5, 1.0};
  Matrix moment = (1.0/3.0) * init_matrix<K>(2,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each segment, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Segment& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT delta[4] = {t[0].x(), t[1].x(), 
		   t[0].y(), t[1].y()};
    Matrix transformation = init_matrix<K>(2,delta);
    FT length = std::sqrt(t.squared_length());
    CGAL_assertion(length != 0.0);

    // Find the 2nd order moment for the segment wrt to the origin by an affine transformation.
    
    // Transform the standard 2nd order moment using the transformation matrix
    transformation = length * transformation * moment * LA::transpose(transformation);
    
    // add to covariance matrix
    covariance[0] += transformation[0][0];
    covariance[1] += transformation[0][1];
    covariance[2] += transformation[1][1];

    mass += length;
  }

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += mass * (-1.0 * c.x() * c.x());
  covariance[1] += mass * (-1.0 * c.x() * c.y());
  covariance[2] += mass * (-1.0 * c.y() * c.y());

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in ascending order, 
  // eigen vectors are sorted in accordance.
  typename DiagonalizeTraits::Vector eigen_values = {{ 0. , 0. }};
  typename DiagonalizeTraits::Matrix eigen_vectors = {{ 0., 0., 0. }};
  DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
    (covariance, eigen_values, eigen_vectors);

  // check unicity and build fitting line accordingly
  if(eigen_values[0] != eigen_values[1])
  {
    // regular case
    line = Line(c, Vector(eigen_vectors[2],eigen_vectors[3]));
    return (FT)1.0 - eigen_values[0] / eigen_values[1];
  } 
  else
  {
    // isotropic case (infinite number of directions)
    // by default: assemble a line that goes through 
    // the centroid and with a default horizontal vector.
    line = Line(c, Vector(1.0, 0.0));
    return (FT)0.0;
  } 
} // end linear_least_squares_fitting_2 for segment set with 1D tag

template < typename InputIterator, typename K, typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Segment_2*,// used for indirection
                               const K& k,                   // kernel
			       const CGAL::Dimension_tag<0>& tag,
			       const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename K::Point_2  Point;
  typedef typename K::Segment_2 Segment;
 
  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  std::list<Point> points;  
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Segment& s = *it;
    points.push_back(s[0]);
    points.push_back(s[1]);
  } 
  return linear_least_squares_fitting_2(points.begin(),points.end(),line,c,k,(Point*)NULL,tag,
					diagonalize_traits);

} // end linear_least_squares_fitting_2 for segment set with 1D tag

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_2_H
