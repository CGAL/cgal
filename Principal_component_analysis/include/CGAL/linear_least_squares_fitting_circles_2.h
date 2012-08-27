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
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_CIRCLES_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_CIRCLES_2_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/centroid.h>
#include <CGAL/eigen_2.h>
#include <CGAL/eigen.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/PCA_util.h>

#include <iterator>
#include <vector>
#include <cmath>

namespace CGAL {

namespace internal {
// Fits a line to a 2D circle set.
// Returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best  (zero variance orthogonally to the fitting line);
//  0 is worst (isotropic case, returns a line with horizontal
//              direction by default)

template < typename InputIterator, typename K >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Circle_2*,// used for indirection
                               const K&,                   // kernel
			                         const CGAL::Dimension_tag<2>& tag)
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Line_2   Line;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Circle_2 Circle;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // ::::::::::DISK:::::::::::::::
  // compute centroid
  c = centroid(first,beyond,K(),tag);

  // assemble covariance matrix as a semi-definite matrix. 
  // Matrix numbering:
  // 0
  // 1 2
  //Final combined covariance matrix for all circles and their combined mass
  FT mass = 0.0;
  FT covariance[3] = {0.0,0.0,0.0};

  // assemble 2nd order moment about the origin.  
  FT temp[4] = {0.25, 0.0,
		0.0, 0.25};
  Matrix moment = init_matrix<K>(2,temp);
  //  Matrix moment = Matrix(2,true,PI);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each circle, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Circle& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT radius = std::sqrt(t.squared_radius());
    FT delta[4] = {radius, 0.0, 
		   0.0, radius};
    Matrix transformation = init_matrix<K>(2,delta);
    FT area = t.squared_radius();
    CGAL_assertion(area != 0.0);

    // Find the 2nd order moment for the circle wrt to the origin by an affine transformation.
    
    // Transform the standard 2nd order moment using the transformation matrix
    transformation = area * transformation * moment * LA::transpose(transformation);
    
    // Translate the 2nd order moment to the center of the circle.
    FT x0 = t.center().x();
    FT y0 = t.center().y();

    // and add to covariance matrix
    covariance[0] += transformation[0][0] + area * x0*x0;
    covariance[1] += transformation[0][1] + area * x0*y0;
    covariance[2] += transformation[1][1] + area * y0*y0;

    mass += area;
  }

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += mass * (-1.0 * c.x() * c.x());
  covariance[1] += mass * (-1.0 * c.x() * c.y());
  covariance[2] += mass * (-1.0 * c.y() * c.y());

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  std::pair<FT,FT> eigen_values;
  std::pair<Vector,Vector> eigen_vectors;
  //  internal::eigen_symmetric_2<K>(final_cov, eigen_vectors, eigen_values);
  FT eigen_vectors1[4];
  FT eigen_values1[2];
  eigen_symmetric<FT>(covariance,2, eigen_vectors1, eigen_values1);
  eigen_values = std::make_pair(eigen_values1[0],eigen_values1[1]);
  eigen_vectors = std::make_pair(Vector(eigen_vectors1[0],eigen_vectors1[1]),Vector(eigen_vectors1[2],eigen_vectors1[3]));
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
  
} // end linear_least_squares_fitting_2 for circle set with 2D tag

template < typename InputIterator, typename K >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Circle_2*,// used for indirection
                               const K&,                   // kernel
			                         const CGAL::Dimension_tag<1>& tag)
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Line_2   Line;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Circle_2 Circle;
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
  //Final combined covariance matrix for all circles and their combined mass
  FT mass = 0.0;
  FT covariance[3] = {0.0,0.0,0.0};

  // assemble 2nd order moment about the origin.  
  FT temp[4] = {1.0, 0.0,
		0.0, 1.0};
  Matrix moment = init_matrix<K>(2,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each circle, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Circle& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT radius = std::sqrt(t.squared_radius());
    FT delta[4] = {radius, 0.0, 
		   0.0, radius};
    Matrix transformation = init_matrix<K>(2,delta);
    FT length = 2 * radius;
    CGAL_assertion(length != 0.0);

    // Find the 2nd order moment for the circle wrt to the origin by an affine transformation.
    
    // Transform the standard 2nd order moment using the transformation matrix
    transformation = 0.5 * length * transformation * moment * LA::transpose(transformation);
    
    // Translate the 2nd order moment to the center of the circle.
    FT x0 = t.center().x();
    FT y0 = t.center().y();

    // and add to covariance matrix
    covariance[0] += transformation[0][0] + length * x0*x0;
    covariance[1] += transformation[0][1] + length * x0*y0;
    covariance[2] += transformation[1][1] + length * y0*y0;

    mass += length;
  }

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += mass * (-1.0 * c.x() * c.x());
  covariance[1] += mass * (-1.0 * c.x() * c.y());
  covariance[2] += mass * (-1.0 * c.y() * c.y());

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  std::pair<FT,FT> eigen_values;
  std::pair<Vector,Vector> eigen_vectors;
  //  internal::eigen_symmetric_2<K>(final_cov, eigen_vectors, eigen_values);
  FT eigen_vectors1[4];
  FT eigen_values1[2];
  eigen_symmetric<FT>(covariance,2, eigen_vectors1, eigen_values1);
  eigen_values = std::make_pair(eigen_values1[0],eigen_values1[1]);
  eigen_vectors = std::make_pair(Vector(eigen_vectors1[0],eigen_vectors1[1]),Vector(eigen_vectors1[2],eigen_vectors1[3]));
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
} // end linear_least_squares_fitting_2 for circle set with 1D tag

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_CIRCLES_2_H
