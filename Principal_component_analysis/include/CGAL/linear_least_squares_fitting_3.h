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

// complete set of parameters
template < typename InputIterator, 
           typename Object,
           typename Kernel,
           typename Tag >
inline
typename Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename Object& object, // plane or line (result)
                               typename Kernel::Point_3& centroid,
                               const Tag& tag,
			                         const Kernel& kernel)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return CGALi::linear_least_squares_fitting_3(first, beyond, object,
                                               centroid, (Value_type*) NULL, kernel, tag);
}


// deduces kernel from value type of input iterator
template < typename InputIterator, 
           typename Object,
					 typename Point,
           typename Tag>
inline
typename Kernel_traits<Object>::Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object,
                               Point& centroid,
			                         const Tag& tag)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return CGAL::linear_least_squares_fitting_3(first,beyond,object,centroid,tag,Kernel());
}

// does not writes centroid and deduces kernel
template < typename InputIterator, 
           typename Object,
           typename Tag>
inline
typename Kernel_traits<Object>::Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object,
			                         const Tag& tag)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  typename Kernel::Point_3 centroid; // unused by caller
  return CGAL::linear_least_squares_fitting_3(first,beyond,object,centroid,tag);
}

// does not return the centroid, deduces kernel, and default tag
// simplest (default) call of the function
/*
template < typename InputIterator, 
           typename Object>
inline
typename Kernel_traits<Object>::Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  typename Kernel::Point_3 centroid; // unused by caller
  return CGAL::linear_least_squares_fitting_3(first,beyond,object,centroid,(Value_type *)NULL,K());
}*/

/*
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

	// default tag (dimension of fitted objects, ie 1 for segments, 2 for triangles, 3 for tets, etc.)
template < typename InputIterator, 
           typename Object,
           typename Kernel >
inline
typename Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename Object& object, // plane or line
                               typename Kernel::Point_3& centroid,
			                         const Kernel& kernel)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return CGALi::linear_least_squares_fitting_3(first, beyond, object,
                                               centroid, (Value_type*) NULL, kernel);
}



// omits centroid
template < typename InputIterator, 
           typename Object,
           typename Kernel,
           typename Tag >
inline
typename Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object,
			                         const Kernel& kernel,
                               const Tag& tag)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typename Kernel::Point_3 centroid; // unused by caller
  return CGALi::linear_least_squares_fitting_3(first, beyond, object,
                                               centroid, (Value_type*) NULL, kernel, tag);
}


	
	*/


CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_3_H
