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

// Functions to compute the centroid of N points.
// Works in 2D and 3D.

// TODO : Note : more numerically stable variants could be implemented as well.
// TODO : Specify a traits class concept ?
// TODO : Grep for "barycenter" and "centroid" in CGAL to check existing usages.
// TODO : Add barycentric_coordinates() (to the kernel, this time).

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// computes the centroid of a 2D point set
// takes an iterator range over K::Point_2
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K&,
         const typename K::Point_2*)
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

// computes the centroid of a 3D point set
// takes an iterator range over K::Point_3
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Point_3*)
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
}// end centroid of a 3D point set

// computes the centroid of a 2D segment set
// takes an iterator range over K::Segment_2
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Segment_2*)
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
} // end centroid of a 2D segment set

// computes the centroid of a 2D triangle set
// takes an iterator range over K::Triangle_2
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Triangle_2*)
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
} // end centroid of a 2D triangle set

// computes the centroid of a 3D triangle set
// takes an iterator range over K::Triangle_3
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& ,
         const typename K::Triangle_3*)
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
} // end centroid of a 3D triangle set

} // namespace CGALi

// computes the centroid of a set of kernel objects
// takes an iterator range over kernel objects
template < typename InputIterator, 
           typename K >
inline
typename Point<Dimension<typename std::iterator_traits<InputIterator>::value_type, K>::value,
               K
              >::type
centroid(InputIterator begin,
         InputIterator end, 
         const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return CGALi::centroid(begin, end, k,(Value_type*) NULL);
}

// this one takes an iterator range over kernel objects
// and uses Kernel_traits<> to find out its kernel.
template < typename InputIterator >
inline
typename Point<Dimension<
               typename std::iterator_traits<InputIterator>::value_type,
               typename Kernel_traits<typename std::iterator_traits<InputIterator>::value_type>::Kernel >::value,
               typename Kernel_traits<typename std::iterator_traits<InputIterator>::value_type>::Kernel >::type
centroid(InputIterator begin, InputIterator end)
{
  typedef typename std::iterator_traits<InputIterator>::value_type  Point;
  typedef typename Kernel_traits<Point>::Kernel                     K;
  return CGAL::centroid(begin, end, K());
}

CGAL_END_NAMESPACE

#endif
