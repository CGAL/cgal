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
// Author(s)     : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_3_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_3_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/IO/io.h>
#include <CGAL/linear_least_squares_fitting_points_3.h>
#include <CGAL/linear_least_squares_fitting_segments_3.h>
#include <CGAL/linear_least_squares_fitting_triangles_3.h>
#include <CGAL/linear_least_squares_fitting_spheres_3.h>
#include <CGAL/linear_least_squares_fitting_cuboids_3.h>
#include <CGAL/linear_least_squares_fitting_tetrahedra_3.h>
#include <CGAL/PCA_tags.h>

#include <iterator>
#include <list>
#include <string>

CGAL_BEGIN_NAMESPACE

template < typename InputIterator, 
           typename K , typename tag>
inline
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,
                               typename K::Point_3& centroid,
                               const K& k,
			                         const tag& t)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return CGALi::linear_least_squares_fitting_3(first, beyond, plane,
                                               centroid, k, (Value_type*) NULL, t);
}

template < typename InputIterator, 
           typename K, typename tag >
inline
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,
                               typename K::Point_3& centroid,
                               const K& k,
			                         const tag& t)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  //  BOOST_STATIC_ASSERT((boost::is_same<typename CGAL::Algebraic_structure_traits<Value_type>::Algebraic_category,CGAL::Field_with_sqrt_tag>::value));
  return CGALi::linear_least_squares_fitting_3(first, beyond, line,
                                               centroid, k, (Value_type*) NULL, t);
}

template < typename InputIterator, 
           typename K, typename tag >
inline
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,
                               const K& k,
			                         const tag& t)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  //   BOOST_STATIC_ASSERT((boost::is_same<typename CGAL::Algebraic_structure_traits<Value_type>::Algebraic_category,CGAL::Field_with_sqrt_tag>::value));
  typename K::Point_3 centroid;
  return CGALi::linear_least_squares_fitting_3(first, beyond, plane,
                                               centroid, k,(Value_type*) NULL, t);
}

template < typename InputIterator, 
           typename K, typename tag >
inline
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,
                               const K& k,
			                         const tag& t)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typename K::Point_3 centroid;
  return CGALi::linear_least_squares_fitting_3(first, beyond, line,
                                               centroid, k,(Value_type*) NULL, t);
}

// deduces kernel 
template < typename InputIterator, 
           typename Object,
           typename tag>
inline
typename Kernel_traits<Object>::Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object,
                               typename Kernel_traits<Object>::Kernel::Point_3& centroid,
			                         const tag& t)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel K;
  return CGAL::linear_least_squares_fitting_3(first,beyond,object,centroid,K(), t);
}

/*
// does not return the centroid and deduces the kernel as well.
template < typename InputIterator, 
           typename Object, typename tag >
inline
typename Kernel_traits<Object>::Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object,
			                         const tag& t)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  //   BOOST_STATIC_ASSERT((boost::is_same<typename CGAL::Algebraic_structure_traits<Value_type>::Algebraic_category,CGAL::Field_with_sqrt_tag>::value));
  typedef typename Kernel_traits<Value_type>::Kernel K;
  return CGAL::linear_least_squares_fitting_3(first,beyond,object,K(), t);
}

// default tag
template < typename InputIterator, 
           typename Object >
inline
typename Kernel_traits<Object>::Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel K;
  return CGAL::linear_least_squares_fitting_3(first,beyond,object,K());
}
*/

// does not return the centroid and deduces the kernel as well.
template < typename InputIterator, 
           typename Object, 
           typename tag = typename PCA_default_dimension< typename std::iterator_traits<InputIterator>::value_type >::Tag>
inline
typename Kernel_traits<Object>::Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object,
			                         const tag& t = tag())
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel K;
  return CGAL::linear_least_squares_fitting_3(first,beyond,object,K(), t);
}

  //   BOOST_STATIC_ASSERT((boost::is_same<typename CGAL::Algebraic_structure_traits<Value_type>::Algebraic_category,CGAL::Field_with_sqrt_tag>::value));


CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_3_H
