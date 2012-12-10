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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_3_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_3_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/centroid.h>
#include <CGAL/eigen.h>
#include <CGAL/PCA_util.h>
#include <CGAL/linear_least_squares_fitting_points_3.h>

#include <iterator>

namespace CGAL {

namespace internal {

// fits a plane to a 3D segment set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Segment_3*,  // used for indirection
                               const K& k,                   // kernel
			                         const CGAL::Dimension_tag<1>& tag)
{
  typedef typename K::FT          FT;
  typedef typename K::Segment_3  Segment;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  // compute centroid
  c = centroid(first,beyond,k,tag);

  // assemble covariance matrix
  FT covariance[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Segment*) NULL,tag);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k);

} // end linear_least_squares_fitting_segments_3

// fits a plane to a 3D segment set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Segment_3*, // used for indirection
                               const K& k,                   // kernel
			                         const CGAL::Dimension_tag<0>& tag)
{
  typedef typename K::Segment_3  Segment;
  typedef typename K::Point_3  Point;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  std::list<Point> points;
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Segment& t = *it;
    points.push_back(t[0]);
    points.push_back(t[1]);
  }

  // compute fitting plane
  return linear_least_squares_fitting_3(points.begin(),points.end(),plane,c,(Point*)NULL,k,tag);

} // end linear_least_squares_fitting_segments_3

// fits a line to a 3D segment set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,      // best fit line
                               typename K::Point_3& c,        // centroid
                               const typename K::Segment_3*,  // used for indirection
                               const K& k,                    // kernel
			                         const CGAL::Dimension_tag<1>& tag)
{
  typedef typename K::FT          FT;
  typedef typename K::Segment_3  Segment;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  // compute centroid
  c = centroid(first,beyond,k,tag);

  // assemble covariance matrix
  FT covariance[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Segment*) NULL,tag);

  // compute fitting line
  return fitting_line_3(covariance,c,line,k);
  
} // end linear_least_squares_fitting_segments_3

// fits a plane to a 3D segment set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,      // best fit line
                               typename K::Point_3& c,        // centroid
                               const typename K::Segment_3*,  // used for indirection
                               const K& k,                    // kernel
			                         const CGAL::Dimension_tag<0>& tag)
{
  typedef typename K::Segment_3  Segment;
  typedef typename K::Point_3  Point;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  std::list<Point> points;
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Segment& t = *it;
    points.push_back(t[0]);
    points.push_back(t[1]);
  }

  // compute fitting plane
  return linear_least_squares_fitting_3(points.begin(),points.end(),line,c,(Point*)NULL,k,tag);

} // end linear_least_squares_fitting_segments_3

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_3_H
