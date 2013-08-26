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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H

#include <CGAL/basic.h>
//#include <CGAL/Algebraic_structure_traits.h>
//#include <CGAL/IO/io.h>
#include <CGAL/linear_least_squares_fitting_points_2.h>
#include <CGAL/linear_least_squares_fitting_segments_2.h>
#include <CGAL/linear_least_squares_fitting_triangles_2.h>
#include <CGAL/linear_least_squares_fitting_circles_2.h>
#include <CGAL/linear_least_squares_fitting_rectangles_2.h>
#include <CGAL/Dimension.h>

#include <iterator>

namespace CGAL {

template < typename InputIterator, 
           typename Kernel,
					 typename Tag>
inline
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename Kernel::Line_2& line,
                               typename Kernel::Point_2& centroid,
                               const Tag& tag,
															 const Kernel& kernel)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return internal::linear_least_squares_fitting_2(first, beyond, line,
                                               centroid,(Value_type*)NULL,kernel,tag);
}

// deduces the kernel from the points in container.
template < typename InputIterator, 
           typename Line,
           typename Point,
	         typename Tag>
inline
typename Kernel_traits<Line>::Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               Line& line,
                               Point& centroid,
			                         const Tag& tag)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
	return CGAL::linear_least_squares_fitting_2(first,beyond,line,centroid,tag,Kernel());
}


template < typename InputIterator, 
           typename Line,
					 typename Tag >
inline
typename Kernel_traits<Line>::Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               Line& line,
                               const Tag& tag)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  typename Kernel::Point_2 centroid; // unused
  return CGAL::linear_least_squares_fitting_2(first,beyond,line,centroid,tag,Kernel());
}


} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H
