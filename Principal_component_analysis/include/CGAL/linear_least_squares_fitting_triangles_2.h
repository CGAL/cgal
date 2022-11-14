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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_2_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
#include <CGAL/centroid.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/PCA_util.h>

#include <iterator>
#include <vector>
#include <cmath>

namespace CGAL {

namespace internal {
// Fits a line to a 2D triangle set.
// Returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best  (zero variance orthogonally to the fitting line);
//  0 is worst (isotropic case, returns a line with horizontal
//              direction by default)

template < typename InputIterator,
           typename Kernel, typename DiagonalizeTraits >
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename Kernel::Line_2& line,   // best fit line
                               typename Kernel::Point_2& c,     // centroid
                               const typename Kernel::Triangle_2*,// used for indirection
                               const Kernel&,                   // kernel
                               const CGAL::Dimension_tag<2>& tag,
                               const DiagonalizeTraits&)
{
  // types
  typedef typename Kernel::FT       FT;
  typedef typename Kernel::Line_2   Line;
  typedef typename Kernel::Vector_2 Vector;
  typedef typename Kernel::Triangle_2 Triangle;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,Kernel(),tag);

  // assemble covariance matrix as a semi-definite matrix.
  // Matrix numbering:
  // 0
  // 1 2
  //Final combined covariance matrix for all triangles and their combined mass
  FT mass = 0.0;
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0. }};

  // assemble the 2nd order moment about the origin.
  FT temp[4] = {FT(1/12.0), FT(1/24.0),
                FT(1/24.0), FT(1/12.0)};

  Matrix moment = init_matrix<FT>(2,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each triangle, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Triangle& t = *it;

    // defined for convenience.
    FT x0 = t[0].x();
    FT y0 = t[0].y();
    FT delta[4] = {t[1].x() - x0, t[2].x() - x0,
                   t[1].y() - y0, t[2].y() - y0};

    Matrix transformation = init_matrix<FT>(2,delta);
    FT area = FT(0.5) * CGAL::abs(LA::determinant(transformation));
    CGAL_assertion(!CGAL::is_zero(area));

    // Find the 2nd order moment for the triangle wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = 2 * area * transformation * moment * LA::transpose(transformation);

    // Translate the 2nd order moment to (x0,y0).
    FT xav0 = (delta[0]+delta[1])/FT(3);
    FT yav0 = (delta[2]+delta[3])/FT(3);

    // and add to the covariance matrix
    covariance[0] += transformation[0][0] + area * (x0*xav0*2 + CGAL::square(x0));
    covariance[1] += transformation[0][1] + area * (x0*yav0 + xav0*y0 + x0*y0);
    covariance[2] += transformation[1][1] + area * (y0*yav0*2 + CGAL::square(y0));

    mass += area;
  }

  CGAL_assertion_msg (!CGAL::is_zero(mass), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] -= mass * (CGAL::square(c.x()));
  covariance[1] -= mass * (c.x() * c.y());
  covariance[2] -= mass * (CGAL::square(c.y()));

  //  std::cout<<"cov: "<<covariance[0]*covariance[2]<<" =? "<<covariance[1]*covariance[1]<<std::endl;

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
} // end linear_least_squares_fitting_2 for triangle set with 2D tag

template < typename InputIterator,
           typename Kernel,
           typename DiagonalizeTraits >
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename Kernel::Line_2& line,   // best fit line
                               typename Kernel::Point_2& c,     // centroid
                               const typename Kernel::Triangle_2*,// used for indirection
                               const Kernel&,                   // kernel
                               const CGAL::Dimension_tag<1>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename Kernel::Triangle_2 Triangle;
  typedef typename Kernel::Segment_2  Segment;
  auto converter = [](const Triangle& t, int idx) -> Segment { return Segment(t[idx], t[(idx+1)%3]); };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_2
    (make_subiterator<Segment, 3> (first, converter),
     make_subiterator<Segment, 3> (beyond),
     line,c,(Segment*)nullptr,Kernel(),tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_2 for triangle set with 1D tag

template < typename InputIterator,
           typename Kernel,
           typename DiagonalizeTraits >
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename Kernel::Line_2& line,   // best fit line
                               typename Kernel::Point_2& c,     // centroid
                               const typename Kernel::Triangle_2*,// used for indirection
                               const Kernel&,                   // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types

  typedef typename Kernel::Triangle_2 Triangle;
  typedef typename Kernel::Point_2 Point;
  auto converter = [](const Triangle& t, int idx) -> Point { return t[idx]; };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_2
    (make_subiterator<Point, 3> (first, converter),
     make_subiterator<Point, 3> (beyond),
     line,c,(Point*)nullptr,Kernel(),tag,
     diagonalize_traits);
} // end linear_least_squares_fitting_2 for triangle set with 0D tag

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_2_H
