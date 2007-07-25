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
// $URL$
// $Id$
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/IO/io.h>
#include <CGAL/linear_least_squares_fitting_points_2.h>
#include <CGAL/linear_least_squares_fitting_segments_2.h>
#include <CGAL/linear_least_squares_fitting_triangles_2.h>
#include <CGAL/linear_least_squares_fitting_circles_2.h>
#include <CGAL/linear_least_squares_fitting_rectangles_2.h>
#include <CGAL/PCA_tags.h>

#include <iterator>

CGAL_BEGIN_NAMESPACE

template < typename InputIterator, 
           typename K , typename tag>
inline
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_2& line,
                               typename K::Point_2& centroid,
                               const K& k,
			       const tag& t)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  //  BOOST_STATIC_ASSERT((boost::is_same<typename CGAL::Algebraic_structure_traits<Value_type>::Algebraic_category,CGAL::Field_with_sqrt_tag>::value));
  return CGALi::linear_least_squares_fitting_2(first, beyond, line,
                                               centroid, k, (Value_type*) NULL, t);
}

template < typename InputIterator, 
           typename K, typename tag >
inline
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_2& line,
                               const K& k,
			       const tag& t)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  //  BOOST_STATIC_ASSERT((boost::is_same<typename CGAL::Algebraic_structure_traits<Value_type>::Algebraic_category,CGAL::Field_with_sqrt_tag>::value));
  typename K::Point_2 centroid;
  return CGALi::linear_least_squares_fitting_2(first, beyond, line,
                                               centroid, k,(Value_type*) NULL, t);
}

// deduces the kernel from the points in container.
template < typename InputIterator, 
           typename Line,
	   typename tag>
inline
typename Kernel_traits<Line>::Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               Line& line,
                               typename Kernel_traits<Line>::Kernel::Point_2& centroid,
			       const tag& t)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  //  BOOST_STATIC_ASSERT((boost::is_same<typename CGAL::Algebraic_structure_traits<Value_type>::Algebraic_category,CGAL::Field_with_sqrt_tag>::value));
  typedef typename Kernel_traits<Value_type>::Kernel K;
  return CGAL::linear_least_squares_fitting_2(first,beyond,line,centroid,K(), t);
}

// does not return the centroid and deduces the kernel as well.
template < typename InputIterator, 
           typename Line, typename tag >
inline
typename Kernel_traits<Line>::Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               Line& line,
			       const tag& t)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel K;
  //  BOOST_STATIC_ASSERT((boost::is_same<typename CGAL::Algebraic_structure_traits<Value_type>::Algebraic_category,CGAL::Field_with_sqrt_tag>::value));
  return CGAL::linear_least_squares_fitting_2(first,beyond,line,K(), t);
}

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H
