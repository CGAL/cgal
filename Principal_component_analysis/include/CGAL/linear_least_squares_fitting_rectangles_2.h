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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_RECTANGLES_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_RECTANGLES_2_H

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
// Fits a line to a 2D rectangle set.
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
                               const typename K::Iso_rectangle_2*,// used for indirection
                               const K&,                   // kernel
                               const CGAL::Dimension_tag<2>& tag,
                               const DiagonalizeTraits&)
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Line_2   Line;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Iso_rectangle_2 Iso_rectangle;
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
  //Final combined covariance matrix for all rectangles and their combined mass
  FT mass = 0.0;
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0. }};

  // assemble 2nd order moment about the origin.
  FT temp[4] = {FT(1/3.0), FT(0.25),
                FT(0.25),  FT(1/3.0)};
  Matrix moment = init_matrix<FT>(2,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each rectangle, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Iso_rectangle& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT x0 = t.xmin();
    FT y0 = t.ymin();
    FT x1 = t.xmax();
    FT y2 = t.ymax();

    FT delta[4] = {x1-x0, 0.0,
                   0.0, y2-y0};

    Matrix transformation = init_matrix<FT>(2,delta);
    FT area = (x1-x0)*(y2-y0);

    CGAL_assertion(area != 0.0);

    // Find the 2nd order moment for the rectangle wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = area * transformation * moment * LA::transpose(transformation);

    // Translate the 2nd order moment to the center of the rectangle.
    FT xav0 = (x1-x0)/FT(2);
    FT yav0 = (y2-y0)/FT(2);
    // and add to covariance matrix
    covariance[0] += transformation[0][0] + area * (x0*xav0*2 + x0*x0);
    covariance[1] += transformation[0][1] + area * (x0*yav0 + xav0*y0 + x0*y0);
    covariance[2] += transformation[1][1] + area * (y0*yav0*2 + y0*y0);

    mass += area;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += -mass * (c.x() * c.x());
  covariance[1] += -mass * (c.x() * c.y());
  covariance[2] += -mass * (c.y() * c.y());

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
} // end linear_least_squares_fitting_2 for rectangle set with 2D tag

template < typename InputIterator, typename K, typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Iso_rectangle_2*,// used for indirection
                               const K&,                   // kernel
                               const CGAL::Dimension_tag<1>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename K::Iso_rectangle_2 Iso_rectangle;
  typedef typename K::Segment_2         Segment_2;
  auto converter = [](const Iso_rectangle& r, int idx) -> Segment_2 { return Segment_2(r[idx], r[(idx+1)%4]); };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_2
    (make_subiterator<Segment_2, 4> (first, converter),
     make_subiterator<Segment_2, 4> (beyond),
     line,c,(Segment_2*)nullptr,K(),tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_2 for rectangle set with 1D tag


template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Iso_rectangle_2*,// used for indirection
                               const K&,                   // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename K::Iso_rectangle_2 Iso_rectangle;
  typedef typename K::Point_2         Point_2;
  auto converter = [](const Iso_rectangle& r, int idx) -> Point_2 { return r[idx]; };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_2
    (make_subiterator<Point_2, 4> (first, converter),
     make_subiterator<Point_2, 4> (beyond),
     line,c,(Point_2*)nullptr,K(),tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_2 for rectangle set with 0D tag

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_RECTANGLES_2_H
