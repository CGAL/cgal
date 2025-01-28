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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_3_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_3_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
#include <CGAL/centroid.h>
#include <CGAL/PCA_util.h>
#include <CGAL/Subiterator.h>

#include <list>
#include <iterator>

namespace CGAL {

namespace internal {

// fits a plane to a 3D triangle set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Triangle_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<2>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Triangle_3  Triangle;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K(),tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Triangle*) nullptr,tag, diagonalize_traits);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k,diagonalize_traits);

} // end linear_least_squares_fitting_triangles_3

// fits a plane to a 3D triangle set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Triangle_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<1>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Triangle_3  Triangle;
  typedef typename K::Segment_3  Segment;
  auto converter = [](const Triangle& t, int idx) -> Segment { return Segment(t[idx], t[(idx+1)%3]); };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_3
    (make_subiterator<Segment, 3> (first, converter),
     make_subiterator<Segment, 3> (beyond),
     plane,c,(Segment*)nullptr,k,tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_triangles_3

// fits a plane to a 3D triangle set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Triangle_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Triangle_3  Triangle;
  typedef typename K::Point_3  Point;
  auto converter = [](const Triangle& t, int idx) -> Point { return t[idx]; };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_3
    (make_subiterator<Point, 3> (first, converter),
     make_subiterator<Point, 3> (beyond),
     plane,c,(Point*)nullptr,k,tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_triangles_3

// fits a line to a 3D triangle set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_3& line,     // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Triangle_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<2>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Triangle_3  Triangle;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K(),tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Triangle*) nullptr,tag, diagonalize_traits);

  // compute fitting line
  return fitting_line_3(covariance,c,line,k,diagonalize_traits);

} // end linear_least_squares_fitting_triangles_3

// fits a line to a 3D triangle set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_3& line,   // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Triangle_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<1>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Triangle_3  Triangle;
  typedef typename K::Segment_3  Segment;
  auto converter = [](const Triangle& t, int idx) -> Segment { return Segment(t[idx], t[(idx+1)%3]); };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_3
    (make_subiterator<Segment, 3> (first, converter),
     make_subiterator<Segment, 3> (beyond),
     line,c,(Segment*)nullptr,k,tag,
     diagonalize_traits);
} // end linear_least_squares_fitting_triangles_3

// fits a line to a 3D triangle set
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_3& line,   // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Triangle_3*,  // used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Triangle_3  Triangle;
  typedef typename K::Point_3  Point;
  auto converter = [](const Triangle& t, int idx) -> Point { return t[idx]; };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_3
    (make_subiterator<Point, 3> (first, converter),
     make_subiterator<Point, 3> (beyond),
     line,c,(Point*)nullptr,k,tag,
     diagonalize_traits);
} // end linear_least_squares_fitting_triangles_3

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_3_H
