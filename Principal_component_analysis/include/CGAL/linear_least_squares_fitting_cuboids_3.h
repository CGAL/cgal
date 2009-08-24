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
// $URL: svn+ssh://gankit@scm.gforge.inria.fr/svn/cgal/trunk/Principal_component_analysis/include/CGAL/linear_least_squares_fitting_cuboids_3.h $
// $Id: linear_least_squares_fitting_2.h 37882 2007-04-03 15:15:30Z spion $
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_CUBOIDS_3_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_CUBOIDS_3_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/centroid.h>
#include <CGAL/eigen.h>
#include <CGAL/PCA_util.h>

#include <CGAL/linear_least_squares_fitting_points_3.h>
#include <CGAL/linear_least_squares_fitting_segments_3.h>

#include <iterator>

CGAL_BEGIN_NAMESPACE

namespace internal {

// fits a plane to a 3D cuboid set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
			                         const CGAL::Dimension_tag<3>& tag)
{
  typedef typename K::FT FT;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  // compute centroid
  c = centroid(first,beyond,k,tag);

  // assemble covariance matrix
  FT covariance[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Iso_cuboid*) NULL,tag);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k);


} // end linear_least_squares_fitting_cuboids_3

// fits a plane to a 3D cuboid set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
			                         const CGAL::Dimension_tag<2>& tag)
{
  typedef typename K::FT FT;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  // compute centroid
  c = centroid(first,beyond,k,tag);

  // assemble covariance matrix
  FT covariance[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Iso_cuboid*) NULL,tag);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k);


} // end linear_least_squares_fitting_cuboids_3

// fits a plane to a 3D cuboid set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
			                         const CGAL::Dimension_tag<1>& tag)
{
  typedef typename K::FT FT;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  std::list<Segment> segments;
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Iso_cuboid& t = *it;
    segments.push_back(Segment(t[0],t[1]));
    segments.push_back(Segment(t[1],t[2]));
    segments.push_back(Segment(t[2],t[3]));    
    segments.push_back(Segment(t[3],t[0]));
    segments.push_back(Segment(t[4],t[5]));
    segments.push_back(Segment(t[5],t[6]));
    segments.push_back(Segment(t[6],t[7]));
    segments.push_back(Segment(t[7],t[4]));
    segments.push_back(Segment(t[5],t[0]));
    segments.push_back(Segment(t[1],t[6]));
    segments.push_back(Segment(t[2],t[7]));
    segments.push_back(Segment(t[3],t[4]));
  }

  // compute fitting plane
  return linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,c,(Segment*)NULL,k,tag);

} // end linear_least_squares_fitting_cuboids_3

// fits a plane to a 3D cuboid set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
			                         const CGAL::Dimension_tag<0>& tag)
{
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  std::list<Point> points;
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Iso_cuboid& t = *it;
    points.push_back(t[0]);
    points.push_back(t[1]);
    points.push_back(t[2]);    
    points.push_back(t[3]);
    points.push_back(t[4]);
    points.push_back(t[5]);
    points.push_back(t[6]);    
    points.push_back(t[7]);
  }

  // compute fitting plane
  return linear_least_squares_fitting_3(points.begin(),points.end(),plane,c,(Point*)NULL,k,tag);

} // end linear_least_squares_fitting_cuboids_3

// fits a line to a 3D cuboid set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,     // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
			                         const CGAL::Dimension_tag<3>& tag)
{
  typedef typename K::FT FT;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  // compute centroid
  c = centroid(first,beyond,k,tag);

  // assemble covariance matrix
  FT covariance[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Iso_cuboid*) NULL,tag);
  
  // compute fitting line
  return fitting_line_3(covariance,c,line,k);
  
} // end linear_least_squares_fitting_cuboids_3

// fits a line to a 3D cuboid set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,   // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
			                         const CGAL::Dimension_tag<2>& tag)
{
  typedef typename K::FT FT;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  // compute centroid
  c = centroid(first,beyond,k,tag);

  // assemble covariance matrix
  FT covariance[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Iso_cuboid*) NULL,tag);
  
  // compute fitting line
  return fitting_line_3(covariance,c,line,k);


} // end linear_least_squares_fitting_cuboids_3

// fits a line to a 3D cuboid set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,   // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
			                         const CGAL::Dimension_tag<1>& tag)
{
  typedef typename K::FT FT;
  typedef typename K::Segment_3 Segment;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  std::list<Segment> segments;
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Iso_cuboid& t = *it;
    segments.push_back(Segment(t[0],t[1]));
    segments.push_back(Segment(t[1],t[2]));
    segments.push_back(Segment(t[2],t[3]));    
    segments.push_back(Segment(t[3],t[0]));
    segments.push_back(Segment(t[4],t[5]));
    segments.push_back(Segment(t[5],t[6]));
    segments.push_back(Segment(t[6],t[7]));
    segments.push_back(Segment(t[7],t[4]));
    segments.push_back(Segment(t[5],t[0]));
    segments.push_back(Segment(t[1],t[6]));
    segments.push_back(Segment(t[2],t[7]));
    segments.push_back(Segment(t[3],t[4]));
  }

  // compute fitting line
  return linear_least_squares_fitting_3(segments.begin(),segments.end(),line,c,(Segment*)NULL,k,tag);

} // end linear_least_squares_fitting_cuboids_3

// fits a line to a 3D cuboid set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,   // best fit line
                               typename K::Point_3& c,       // centroid
                               const typename K::Iso_cuboid_3*,  // used for indirection
                               const K& k,                   // kernel
			                         const CGAL::Dimension_tag<0>& tag)
{
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  std::list<Point> points;
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Iso_cuboid& t = *it;
    points.push_back(t[0]);
    points.push_back(t[1]);
    points.push_back(t[2]);    
    points.push_back(t[3]);
    points.push_back(t[4]);
    points.push_back(t[5]);
    points.push_back(t[6]);    
    points.push_back(t[7]);
  }

  // compute fitting line
  return linear_least_squares_fitting_3(points.begin(),points.end(),line,c,(Point*)NULL,k,tag);

} // end linear_least_squares_fitting_cuboids_3

} // end namespace internal

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_CUBOIDS_3_H
