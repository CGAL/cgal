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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_3_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_3_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>

#include <CGAL/linear_least_squares_fitting_points_3.h>
#include <CGAL/linear_least_squares_fitting_segments_3.h>
#include <CGAL/linear_least_squares_fitting_triangles_3.h>
#include <CGAL/linear_least_squares_fitting_cuboids_3.h>
#include <CGAL/linear_least_squares_fitting_tetrahedra_3.h>
#include <CGAL/linear_least_squares_fitting_spheres_3.h>

#include <CGAL/Default_diagonalize_traits.h>

#include <CGAL/Dimension.h>

#include <iterator>
#include <string>

namespace CGAL {

// complete set of parameters
template < typename InputIterator, 
           typename Object,
           typename Kernel,
           typename Tag,
	   typename DiagonalizeTraits >
inline
typename Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object, // plane or line
                               typename Kernel::Point_3& centroid, 
                               const Tag& tag, // dimension tag, ranges from 0 to 3
			       const Kernel& kernel,
			       const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return internal::linear_least_squares_fitting_3(first, beyond, object,
						  centroid, (Value_type*) NULL, kernel, tag,
						  diagonalize_traits);
}

// deduces kernel from value type of input iterator
// use default DiagonalizeTraits
template < typename InputIterator, 
           typename Object,
	   typename Point,
           typename Tag >
inline
typename Kernel_traits<Object>::Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object,  // plane or line
                               Point& centroid,
			       const Tag& tag) // dimension tag, ranges from 0 to 3
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return CGAL::linear_least_squares_fitting_3(first,beyond,object,centroid,tag,Kernel(),
					      Default_diagonalize_traits<typename Kernel::FT, 3>());

}

// deduces kernel and does not write centroid
// use default DiagonalizeTraits
template < typename InputIterator, 
           typename Object,
           typename Tag>
inline
typename Kernel_traits<Object>::Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object, // plane or line
			       const Tag& tag) // dimension tag, ranges from 0 to 3
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  typename Kernel::Point_3 centroid; // not used by caller
  return CGAL::linear_least_squares_fitting_3(first,beyond,object,centroid,tag,Kernel(),
					      Default_diagonalize_traits<typename Kernel::FT, 3>());

}

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_3_H
