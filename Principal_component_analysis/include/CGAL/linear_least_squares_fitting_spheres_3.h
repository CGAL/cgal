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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_SPHERES_3_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_SPHERES_3_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
#include <CGAL/centroid.h>
#include <CGAL/PCA_util.h>

#include <iterator>

namespace CGAL {

namespace internal {

// fits a plane to a set of 3D balls (3D)
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Sphere_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<3>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Sphere_3    Sphere;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K(),tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Sphere*) nullptr,tag, diagonalize_traits);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k,diagonalize_traits);

} // end linear_least_squares_fitting_spheres_3

// fits a plane to a 3D sphere set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Sphere_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<2>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Sphere_3    Sphere;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K(),tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Sphere*) nullptr,tag, diagonalize_traits);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k,diagonalize_traits);

} // end linear_least_squares_fitting_spheres_3


// fits a line to a 3D ball set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_3& line,     // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Sphere_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<3>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Sphere_3  Sphere;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K(),tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Sphere*) nullptr,tag, diagonalize_traits);


  // compute fitting line
  return fitting_line_3(covariance,c,line,k,diagonalize_traits);

} // end linear_least_squares_fitting_spheres_3

// fits a line to a 3D sphere set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_3& line,     // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Sphere_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<2>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Sphere_3  Sphere;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K(),tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Sphere*) nullptr,tag, diagonalize_traits);

  // compute fitting line
  return fitting_line_3(covariance,c,line,k,diagonalize_traits);

} // end linear_least_squares_fitting_spheres_3

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_SPHERES_3_H
