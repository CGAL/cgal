// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_CENTROID_H
#define CGAL_CENTROID_H

#include <CGAL/basic.h>
#include <iterator>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/Dimension_utils.h>
#include <CGAL/Dimension.h>

#include <list>
// Functions to compute the centroid of N points.
// Works in 2D and 3D.

// TODO : Note : more numerically stable variants could be implemented as well.
// TODO : Specify a traits class concept ?
// TODO : Grep for "barycenter" and "centroid" in CGAL to check existing usages.
// TODO : Add barycentric_coordinates() (to the kernel, this time).

CGAL_BEGIN_NAMESPACE

namespace CGALi {

//:::::::::: 2D Objects :::::::::::::::::::

// computes the centroid of a 2D point set
// takes an iterator range over K::Point_2
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K&,
         const typename K::Point_2*,
         CGAL::Dimension_tag<0>)
{
  typedef typename K::Vector_2 Vector;
  typedef typename K::FT FT;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  unsigned int nb_pts = 0;
  while(begin != end) 
  {
    v = v + (*begin++ - ORIGIN);
    nb_pts++;
  }
  return ORIGIN + v / (FT)nb_pts;
}// end centroid of a 2D point set

// computes the centroid of a 2D segment set
// takes an iterator range over K::Segment_2

// centroid for 2D segment set with 0D tag
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Segment_2*,
         CGAL::Dimension_tag<0> tag)
{
  typedef typename K::Point_2  Point;
  typedef typename K::Segment_2 Segment;
  
  CGAL_precondition(begin != end);
  
  std::list<Point> points;  
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Segment& s = *it;
    points.push_back(s[0]);
    points.push_back(s[1]);
  } 
  return centroid(points.begin(),points.end(),k,(Point*)NULL,tag);
}// end centroid for 2D segment set with 0D tag

// centroid for 2D segment set with 1D tag
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Segment_2*,
         CGAL::Dimension_tag<1>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Point_2  Point;
  typedef typename K::Segment_2 Segment;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_lengths = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Segment& s = *it;
    FT length = std::sqrt(std::abs(s.squared_length()));
    Point c = K().construct_midpoint_2_object()(s[0],s[1]);
    v = v + length * (c - ORIGIN);
    sum_lengths += length;
  }
  CGAL_assertion(sum_lengths != 0.0);
  return ORIGIN + v / sum_lengths;
} // end centroid of a 2D segment set with 1D tag

// computes the centroid of a 2D triangle set
// takes an iterator range over K::Triangle_2

// centroid for 2D triangle set with 0D tag
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Triangle_2*,
         CGAL::Dimension_tag<0> tag)
{
  typedef typename K::Triangle_2 Triangle;
  typedef typename K::Point_2 Point;

  CGAL_precondition(begin != end);
  
  std::list<Point> points;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Triangle& triangle = *it;
    points.push_back(triangle[0]);
    points.push_back(triangle[1]);
    points.push_back(triangle[2]);
  }
  return centroid(points.begin(),points.end(),k,(Point*)NULL,tag);

} // end centroid of a 2D triangle set with 0D tag

// centroid for 2D triangle set with 1D tag
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Triangle_2*,
         CGAL::Dimension_tag<1> tag)
{
  typedef typename K::Triangle_2 Triangle;
  typedef typename K::Segment_2 Segment;

  CGAL_precondition(begin != end);
  
  std::list<Segment> segments;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Triangle& triangle = *it;
    segments.push_back(triangle[0],triangle[1]);
    segments.push_back(triangle[1],triangle[2]);
    segments.push_back(triangle[2],triangle[0]);
  }
  return centroid(segments.begin(),segments.end(),k,(Segment*)NULL,tag);

} // end centroid of a 2D triangle set with 1D tag

// centroid for 2D triangle set with 2D tag
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Triangle_2*,
         CGAL::Dimension_tag<2>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Point_2  Point;
  typedef typename K::Triangle_2 Triangle;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_areas = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Triangle& triangle = *it;
    FT unsigned_area = std::abs(triangle.area());
    Point c = K().construct_centroid_2_object()(triangle[0],triangle[1],triangle[2]);
    v = v + unsigned_area * (c - ORIGIN);
    sum_areas += unsigned_area;
  }
  CGAL_assertion(sum_areas != 0.0);
  return ORIGIN + v / sum_areas;
} // end centroid of a 2D triangle set with 2D tag

// computes the centroid of a 2D circle set
// takes an iterator range over K::Circle_2

// centroid for 2D circle set with 1D tag
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Circle_2*,
         CGAL::Dimension_tag<1>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Point_2  Point;
  typedef typename K::Circle_2 Circle;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_lengths = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Circle& s = *it;
    FT radius = std::sqrt(s.squared_radius());
    Point c = s.center();
    v = v + radius * (c - ORIGIN);
    sum_lengths += radius;
  }
  CGAL_assertion(sum_lengths != 0.0);
  return ORIGIN + v / sum_lengths;
} // end centroid of a 2D circle set with 1D tag

// centroid for 2D circle set with 2D tag
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Circle_2*,
         CGAL::Dimension_tag<2>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Point_2  Point;
  typedef typename K::Circle_2 Circle;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_areas = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Circle& s = *it;
    FT sq_radius = s.squared_radius();
    Point c = s.center();
    v = v + sq_radius * (c - ORIGIN);
    sum_areas += sq_radius;
  }
  CGAL_assertion(sum_areas != 0.0);
  return ORIGIN + v / sum_areas;
} // end centroid of a 2D circle set with 2D tag

// computes the centroid of a 2D rectangle set
// takes an iterator range over K::Iso_Rectangle_2

// centroid for 2D rectangle set with 0D tag
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Iso_rectangle_2*,
         CGAL::Dimension_tag<0> tag)
{
  typedef typename K::Iso_rectangle_2 Iso_rectangle;
  typedef typename K::Point_2 Point;

  CGAL_precondition(begin != end);
  
  std::list<Point> points;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Iso_rectangle& r = *it;
    points.push_back(r[0]);
    points.push_back(r[1]);
    points.push_back(r[2]);
    points.push_back(r[3]);
  }
  return centroid(points.begin(),points.end(),k,(Point*)NULL,tag);

} // end centroid of a 2D rectangle set with 0D tag

// centroid for 2D rectangle set with 1D tag
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Iso_rectangle_2*,
         CGAL::Dimension_tag<1> tag)
{
  typedef typename K::Iso_rectangle_2 Iso_rectangle;
  typedef typename K::Segment_2 Segment;

  CGAL_precondition(begin != end);
  
  std::list<Segment> segments;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Iso_rectangle& r = *it;
    segments.push_back(r[0],r[1]);
    segments.push_back(r[1],r[2]);
    segments.push_back(r[2],r[3]);
    segments.push_back(r[3],r[0]);
  }
  return centroid(segments.begin(),segments.end(),k,(Segment*)NULL,tag);

} // end centroid of a 2D rectangle set with 1D tag

// centroid for 2D rectangle set with 2D tag
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Iso_rectangle_2*,
         CGAL::Dimension_tag<2>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Point_2  Point;
  typedef typename K::Iso_rectangle_2 Iso_rectangle;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_areas = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Iso_rectangle& r = *it;
    FT unsigned_area = std::abs(r.area());
    Point c = K().construct_centroid_2_object()(r[0],r[1],r[2],r[3]);
    v = v + unsigned_area * (c - ORIGIN);
    sum_areas += unsigned_area;
  }
  CGAL_assertion(sum_areas != 0.0);
  return ORIGIN + v / sum_areas;
} // end centroid of a 2D rectangle set with 2D tag

// computes the centroid of a 3D point set
// takes an iterator range over K::Point_3

// centroid for 3D point set with 0D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K&,
         const typename K::Point_3*,
         CGAL::Dimension_tag<0>)
{
  typedef typename K::Vector_3 Vector;
  typedef typename K::FT FT;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  unsigned int nb_pts = 0;
  while (begin != end) 
  {
    v = v + (*begin++ - ORIGIN);
    nb_pts++;
  }
  return ORIGIN + v / (FT)nb_pts;
}// end centroid of a 3D point set with 0D tag

// centroid for 3D segment set with 1D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Segment_3*,
         CGAL::Dimension_tag<1>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Point_3  Point;
  typedef typename K::Segment_3 Segment;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_lengths = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Segment& s = *it;
    FT length = std::sqrt(s.squared_length());
		Point c = CGAL::midpoint(s.source(),s.target());
    // Point c = K().construct_midpoint_3_object()(s[0],s[1]);
    //Point c = Point((s[0][0] + s[1][0])/2.0, (s[0][1] + s[1][1])/2.0, (s[0][2] + s[1][2])/2.0);
    v = v + length * (c - ORIGIN);
    sum_lengths += length;
  }
  CGAL_assertion(sum_lengths != 0.0);
  return ORIGIN + v / sum_lengths;
} // end centroid of a 3D segment set with 1D tag

// computes the centroid of a 3D triangle set
// takes an iterator range over K::Triangle_3

// centroid for 3D triangle set with 0D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Triangle_3*,
         CGAL::Dimension_tag<0> tag)
{
  typedef typename K::Triangle_3 Triangle;
  typedef typename K::Point_3 Point;

  CGAL_precondition(begin != end);
  
  std::list<Point> points;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Triangle& triangle = *it;
    points.push_back(triangle[0]);
    points.push_back(triangle[1]);
    points.push_back(triangle[2]);
  }
  return centroid(points.begin(),points.end(),k,(Point*)NULL,tag);

} // end centroid of a 3D triangle set with 0D tag

// centroid for 3D triangle set with 1D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Triangle_3*,
         CGAL::Dimension_tag<1> tag)
{
  typedef typename K::Triangle_3 Triangle;
  typedef typename K::Segment_3 Segment;

  CGAL_precondition(begin != end);
  
  std::list<Segment> segments;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Triangle& triangle = *it;
    segments.push_back(triangle[0],triangle[1]);
    segments.push_back(triangle[1],triangle[2]);
    segments.push_back(triangle[2],triangle[0]);
  }
  return centroid(segments.begin(),segments.end(),k,(Segment*)NULL,tag);

} // end centroid of a 3D triangle set with 1D tag

// centroid for 3D triangle set with 2D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Triangle_3*,
         CGAL::Dimension_tag<2>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Point_3  Point;
  typedef typename K::Triangle_3 Triangle;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_areas = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Triangle& triangle = *it;
    FT unsigned_area = std::sqrt(triangle.squared_area());
    Point c = K().construct_centroid_3_object()(triangle[0],triangle[1],triangle[2]);
    v = v + unsigned_area * (c - ORIGIN);
    sum_areas += unsigned_area;
  }
  CGAL_assertion(sum_areas != 0.0);
  return ORIGIN + v / sum_areas;
} // end centroid of a 3D triangle set with 2D tag

// computes the centroid of a 3D sphere set
// takes an iterator range over K::Sphere_3

// centroid for 3D sphere set with 2D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Sphere_3*,
         CGAL::Dimension_tag<2>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Point_3  Point;
  typedef typename K::Sphere_3 Sphere;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_areas = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Sphere& sphere = *it;
    FT unsigned_area = sphere.squared_radius();
    Point c = sphere.center();
    v = v + unsigned_area * (c - ORIGIN);
    sum_areas += unsigned_area;
  }
  CGAL_assertion(sum_areas != 0.0);
  return ORIGIN + v / sum_areas;
} // end centroid of a 3D sphere set with 2D tag

// centroid for 3D sphere set with 3D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Sphere_3*,
         CGAL::Dimension_tag<3>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Point_3  Point;
  typedef typename K::Sphere_3 Sphere;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_volumes = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Sphere& sphere = *it;
    FT unsigned_volume = sphere.squared_radius() * std::sqrt(sphere.squared_radius());
    Point c = sphere.center();
    v = v + unsigned_volume * (c - ORIGIN);
    sum_volumes += unsigned_volume;
  }
  CGAL_assertion(sum_volumes != 0.0);
  return ORIGIN + v / sum_volumes;
} // end centroid of a 3D sphere set with 3 tag

// computes the centroid of a 3D cuboid set
// takes an iterator range over K::Iso_cuboid_3

// centroid for 3D cuboid set with 0D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Iso_cuboid_3*,
         CGAL::Dimension_tag<0> tag)
{
  typedef typename K::Iso_cuboid_3 Iso_cuboid;
  typedef typename K::Point_3 Point;

  CGAL_precondition(begin != end);
  
  std::list<Point> points;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Iso_cuboid& cuboid = *it;
    points.push_back(cuboid[0]);
    points.push_back(cuboid[1]);
    points.push_back(cuboid[2]);
    points.push_back(cuboid[3]);
    points.push_back(cuboid[4]);
    points.push_back(cuboid[5]);
    points.push_back(cuboid[6]);
    points.push_back(cuboid[7]);
  }
  return centroid(points.begin(),points.end(),k,(Point*)NULL,tag);

} // end centroid of a 3D cuboid set with 0D tag

// centroid for 3D cuboid set with 1D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Iso_cuboid_3*,
         CGAL::Dimension_tag<1> tag)
{
  typedef typename K::Iso_cuboid_3 Iso_cuboid;
  typedef typename K::Segment_3 Segment;

  CGAL_precondition(begin != end);
  
  std::list<Segment> segments;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Iso_cuboid& cuboid = *it;
    segments.push_back(cuboid[0],cuboid[1]);
    segments.push_back(cuboid[1],cuboid[2]);
    segments.push_back(cuboid[2],cuboid[3]);
    segments.push_back(cuboid[3],cuboid[0]);
    segments.push_back(cuboid[0],cuboid[5]);
    segments.push_back(cuboid[5],cuboid[4]);
    segments.push_back(cuboid[4],cuboid[3]);
    segments.push_back(cuboid[1],cuboid[6]);
    segments.push_back(cuboid[6],cuboid[7]);
    segments.push_back(cuboid[7],cuboid[2]);
    segments.push_back(cuboid[4],cuboid[7]);
    segments.push_back(cuboid[5],cuboid[6]);
  }
  return centroid(segments.begin(),segments.end(),k,(Segment*)NULL,tag);

} // end centroid of a 3D cuboid set with 1D tag

// centroid for 3D cuboid set with 2D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Iso_cuboid_3*,
         CGAL::Dimension_tag<2>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Point_3  Point;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_areas = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Iso_cuboid& cuboid = *it;
    FT unsigned_area = 2 * ((cuboid.xmax()-cuboid.xmin())*(cuboid.ymax()-cuboid.ymin()) + (cuboid.xmax()-cuboid.xmin())*(cuboid.zmax()-cuboid.zmin()) + (cuboid.ymax()-cuboid.ymin())*(cuboid.zmax()-cuboid.zmin()));
    Point c = K().construct_centroid_3_object()(cuboid[0],cuboid[1],cuboid[3],cuboid[5]);
    v = v + unsigned_area * (c - ORIGIN);
    sum_areas += unsigned_area;
  }
  CGAL_assertion(sum_areas != 0.0);
  return ORIGIN + v / sum_areas;
} // end centroid of a 3D cuboid set with 2D tag

// centroid for 3D cuboid set with 3D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Iso_cuboid_3*,
         CGAL::Dimension_tag<3>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Point_3  Point;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_volumes = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Iso_cuboid& cuboid = *it;
    FT unsigned_volume = cuboid.volume();
    Point c = K().construct_centroid_3_object()(cuboid[0],cuboid[1],cuboid[3],cuboid[5]);
    v = v + unsigned_volume * (c - ORIGIN);
    sum_volumes += unsigned_volume;
  }
  CGAL_assertion(sum_volumes != 0.0);
  return ORIGIN + v / sum_volumes;
} // end centroid of a 3D cuboid set with 3D tag

// centroid for 3D Tetrahedron set with 3D tag
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Tetrahedron_3*,
         CGAL::Dimension_tag<3>)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Point_3  Point;
  typedef typename K::Tetrahedron_3 Tetrahedron;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_volumes = 0.0;
  for(InputIterator it = begin;
      it != end;
      it++)
  {
    const Tetrahedron& Tetrahedron = *it;
    FT unsigned_volume = Tetrahedron.volume();
    Point c = K().construct_centroid_3_object()(Tetrahedron[0],Tetrahedron[1],Tetrahedron[2],Tetrahedron[3]);
    v = v + unsigned_volume * (c - ORIGIN);
    sum_volumes += unsigned_volume;
  }
  CGAL_assertion(sum_volumes != 0.0);
  return ORIGIN + v / sum_volumes;
} // end centroid of a 3D Tetrahedron set with 3D tag

} // namespace CGALi

// computes the centroid of a set of kernel objects
// takes an iterator range over kernel objects
template < typename InputIterator, 
           typename K, 
           typename Dim_tag  >
inline
typename Point<typename Ambiant_dimension<typename std::iterator_traits<InputIterator>::value_type, K>::type,
               K
              >::type
centroid(InputIterator begin,
         InputIterator end, 
         const K& k,
         Dim_tag tag)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return CGALi::centroid(begin, end, k,(Value_type*) NULL, tag);
}

// this one takes an iterator range over kernel objects
// and uses Kernel_traits<> to find out its kernel.
template < typename InputIterator, typename Dim_tag >
inline
typename Point<typename Ambiant_dimension<
               typename std::iterator_traits<InputIterator>::value_type,
               typename Kernel_traits<typename std::iterator_traits<InputIterator>::value_type>::Kernel >::type,
               typename Kernel_traits<typename std::iterator_traits<InputIterator>::value_type>::Kernel >::type
centroid(InputIterator begin, InputIterator end, Dim_tag tag)
{
  typedef typename std::iterator_traits<InputIterator>::value_type  Point;
  typedef typename Kernel_traits<Point>::Kernel                     K;
  return CGAL::centroid(begin, end, K(), tag);
}

CGAL_END_NAMESPACE

#endif
