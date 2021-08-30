// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
//#include <CGAL/Algebraic_structure_traits.h>
//#include <CGAL/IO/io.h>
#include <CGAL/linear_least_squares_fitting_points_2.h>
#include <CGAL/linear_least_squares_fitting_segments_2.h>
#include <CGAL/linear_least_squares_fitting_triangles_2.h>
#include <CGAL/linear_least_squares_fitting_circles_2.h>
#include <CGAL/linear_least_squares_fitting_rectangles_2.h>
#include <CGAL/Dimension.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <iterator>

namespace CGAL {

template < typename InputIterator,
           typename Kernel,
           typename Tag,
           typename DiagonalizeTraits>
inline
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename Kernel::Line_2& line,
                               typename Kernel::Point_2& centroid,
                               const Tag& tag,
                               const Kernel& kernel,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return internal::linear_least_squares_fitting_2(first, beyond, line,
                                                  centroid,(Value_type*)nullptr,kernel,tag,
                                                  diagonalize_traits);
}

// deduces the kernel from the points in container.
// Use default DiagonalizeTraits
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
  return CGAL::linear_least_squares_fitting_2
    (first,beyond,line,centroid,tag,Kernel(),
     Default_diagonalize_traits<typename Kernel::FT, 2>());
}

template < typename InputIterator,
           typename Line,
           typename Tag>
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
  return CGAL::linear_least_squares_fitting_2(first,beyond,line,centroid,tag,Kernel(),
                                              Default_diagonalize_traits<typename Kernel::FT, 2>());
}


} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_2_H
