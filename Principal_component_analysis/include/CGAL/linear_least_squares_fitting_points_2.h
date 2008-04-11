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
// $URL: svn+ssh://gankit@scm.gforge.inria.fr/svn/cgal/trunk/Principal_component_analysis/include/CGAL/linear_least_squares_fitting_triangles.h $
// $Id: linear_least_squares_fitting_2.h 37882 2007-04-03 15:15:30Z spion $
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_POINTS_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_POINTS_2_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/centroid.h>
#include <CGAL/eigen_2.h>

#include <iterator>
#include <cmath>

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
                               const typename K::Point_2*,// used for indirection
                               const K&,                   // kernel
			                         const CGAL::Dimension_tag<0>& tag)
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Line_2   Line;
  typedef typename K::Point_2  Point;
  typedef typename K::Vector_2 Vector;

	// if internally double, declare a kernel

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K(),tag);

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

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_POINTS_2_H
