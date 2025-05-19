// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_2_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
#include <CGAL/centroid.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/PCA_util.h>
#include <CGAL/Subiterator.h>

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
  Matrix moment = FT(1.0/3.0) * init_matrix<FT>(2,temp);

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
    Matrix transformation = init_matrix<FT>(2,delta);
    using std::sqrt;
    FT length = CGAL::approximate_sqrt(t.squared_length());
    CGAL_assertion(!CGAL::is_zero(length));

    // Find the 2nd order moment for the segment wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = length * transformation * moment * LA::transpose(transformation);

    // add to covariance matrix
    covariance[0] += transformation[0][0];
    covariance[1] += transformation[0][1];
    covariance[2] += transformation[1][1];

    mass += length;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] -= mass * CGAL::square(c.x());
  covariance[1] -= mass * (c.x() * c.y());
  covariance[2] -= mass * (CGAL::square(c.y()));

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
    line = Line(c, Vector(FT(1), FT(0)));
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
  auto converter = [](const Segment& s, int idx) -> Point { return s[idx]; };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_2
    (make_subiterator<Point, 2> (first, converter),
     make_subiterator<Point, 2> (beyond),
     line,c,(Point*)nullptr,k,tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_2 for segment set with 1D tag

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_2_H
