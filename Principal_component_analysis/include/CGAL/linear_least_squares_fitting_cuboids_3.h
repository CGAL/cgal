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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_CUBOIDS_3_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_CUBOIDS_3_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
#include <CGAL/centroid.h>
#include <CGAL/PCA_util.h>
#include <CGAL/linear_least_squares_fitting_points_3.h>
#include <CGAL/linear_least_squares_fitting_segments_3.h>
#include <CGAL/Subiterator.h>

#include <list>
#include <iterator>

namespace CGAL {

namespace internal {

// fits a plane to a 3D cuboid set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<3>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,k,tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Iso_cuboid*) nullptr,tag, diagonalize_traits);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k,diagonalize_traits);


} // end linear_least_squares_fitting_cuboids_3

// fits a plane to a 3D cuboid set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<2>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,k,tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Iso_cuboid*) nullptr,tag,diagonalize_traits);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k,diagonalize_traits);


} // end linear_least_squares_fitting_cuboids_3

// fits a plane to a 3D cuboid set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<1>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Segment_3 Segment;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  auto converter = [](const Iso_cuboid& c, int idx) -> Segment
    {
      if (idx < 7)
        return Segment (c[idx], c[idx+1]);
      if (idx < 10)
        return Segment (c[idx - 6], c[(idx-1)%8]);
      if (idx == 10)
        return Segment (c[5], c[0]);
      CGAL_assertion (idx == 11);
      return Segment (c[7], c[4]);
    };

  return linear_least_squares_fitting_3
    (make_subiterator<Segment, 12> (first, converter),
     make_subiterator<Segment, 12> (beyond),
     plane,c,(Segment*)nullptr,k,tag,diagonalize_traits);

} // end linear_least_squares_fitting_cuboids_3

// fits a plane to a 3D cuboid set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Point_3 Point;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;
  auto converter = [](const Iso_cuboid& c, int idx) -> Point { return c[idx]; };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_3
    (make_subiterator<Point, 8> (first, converter),
     make_subiterator<Point, 8> (beyond),
     plane,c,(Point*)nullptr,k,tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_cuboids_3

// fits a line to a 3D cuboid set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_3& line,     // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<3>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,k,tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Iso_cuboid*) nullptr,tag,diagonalize_traits);

  // compute fitting line
  return fitting_line_3(covariance,c,line,k,diagonalize_traits);

} // end linear_least_squares_fitting_cuboids_3

// fits a line to a 3D cuboid set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_3& line,   // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<2>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,k,tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Iso_cuboid*) nullptr,tag, diagonalize_traits);

  // compute fitting line
  return fitting_line_3(covariance,c,line,k,diagonalize_traits);


} // end linear_least_squares_fitting_cuboids_3

// fits a line to a 3D cuboid set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_3& line,   // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<1>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Segment_3 Segment;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;
  auto converter = [](const Iso_cuboid& c, int idx) -> Segment
    {
      if (idx < 7)
        return Segment (c[idx], c[idx+1]);
      if (idx < 10)
        return Segment (c[idx - 6], c[(idx-1)%8]);
      if (idx == 10)
        return Segment (c[5], c[0]);
      CGAL_assertion (idx == 11);
      return Segment (c[7], c[4]);
    };


  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_3
    (make_subiterator<Segment, 12> (first, converter),
     make_subiterator<Segment, 12> (beyond),
     line,c,(Segment*)nullptr,k,tag,diagonalize_traits);

} // end linear_least_squares_fitting_cuboids_3

// fits a line to a 3D cuboid set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_3& line,   // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Point_3 Point;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;
  auto converter = [](const Iso_cuboid& c, int idx) -> Point { return c[idx]; };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_3
    (make_subiterator<Point, 8> (first, converter),
     make_subiterator<Point, 8> (beyond),
     line,c,(Point*)nullptr,k,tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_cuboids_3

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_CUBOIDS_3_H
