// Copyright (c) 2005  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
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

// computes the centroid of 2D point set
// takes an iterator range over K::Point_2
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Point_2*)
{
  typedef typename K::Vector_2 Vector;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  int nb_pts = 0;
  while (begin != end) 
  {
    v = v + (*begin++ - ORIGIN);
    ++nb_pts;
  }
  return ORIGIN + v / nb_pts;
}

// computes the centroid of 3D point set
// takes an iterator range over K::Point_3
template < typename InputIterator, 
           typename K >
typename K::Point_3
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Point_3*)
{
  typedef typename K::Vector_3 Vector;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  int nb_pts = 0;
  while (begin != end) 
  {
    v = v + (*begin++ - ORIGIN);
    ++nb_pts;
  }
  return ORIGIN + v / nb_pts;
}

// computes the centroid of 2D triangle set
// takes an iterator range over K::Triangle_2
template < typename InputIterator, 
           typename K >
typename K::Point_2
centroid(InputIterator begin, 
         InputIterator end, 
         const K& k,
         const typename K::Triangle_2*)
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Point_2  Point;
  typedef typename K::Triangle_2 Triangle;

  CGAL_precondition(begin != end);

  Vector v = NULL_VECTOR;
  FT sum_area = 0;
  for(InputIterator it = begin;
      it != end;
      it++)
  while (begin != end) 
  {
    const Triangle& triangle = *it;
    FT area = std::abs(triangle.area());
    Point c = centroid(triangle[0],triangle[1],triangle[2]);
    v = v + area * (c - ORIGIN);
    sum_area += area;
  }
  return ORIGIN + v / sum_area;
}

} // namespace CGALi



template < typename InputIterator, 
           typename K >
inline
typename K::Point_2
centroid(InputIterator begin,
         InputIterator end, 
         const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type
    Value_type;
  return CGALi::centroid(begin, end, k,(Value_type*) NULL);
}

// This one takes an iterator range over K::Point_[23]
// And it uses Kernel_traits<> to find out its kernel.
template < typename InputIterator >
inline
typename std::iterator_traits<InputIterator>::value_type
centroid(InputIterator begin, InputIterator end)
{
  typedef typename std::iterator_traits<InputIterator>::value_type  Point;
  typedef typename Kernel_traits<Point>::Kernel                     K;
  return CGAL::centroid(begin, end, K());
}

CGAL_END_NAMESPACE

#endif
