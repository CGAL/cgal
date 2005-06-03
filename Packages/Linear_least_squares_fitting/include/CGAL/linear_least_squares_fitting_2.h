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
// $Source:
// /CVSROOT/CGAL/Packages/Interpolation/include/CGAL/interpolation_functions.h,v $
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Pierre Alliez

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_2
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_2

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <iterator>
#include <list>
#include <string>
#include <CGAL/eigen.h>
#include <CGAL/centroid.h>

CGAL_BEGIN_NAMESPACE

template < typename InputIterator, 
           typename Line, 
           typename K >
void 
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               Line& line,
                               const K& k)
{
  // types
  typedef typename K::FT          FT;
  typedef typename K::Point_2     Point;
  typedef typename K::Vector_2    Vector;
  typedef typename K::Direction_2 Direction;

  // compute centroid
  Point c = centroid(begin,end,K());

  // assemble covariance matrix (semi-definite)
  // numbering:
  // 0          
  // 1 2        
  FT covariance[3] = {0,0,0};
  FT eigen_values[2] = {0,0};
  FT eigen_vectors[4] = {0,0,0,0};
  for(InputIterator it = begin;
		    it != end;
		    it++)
  {
    const Point& p = *it;
    FT dx = p.x()-c.x();
    FT dy = p.y()-c.y();
    covariance[0] += dx * dx;
    covariance[1] += dx * dy;
    covariance[2] += dy * dy;
  }

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  eigen_semi_definite_symmetric(covariance,2,eigen_vectors,eigen_values);

  // check unicity and build fitting line accordingly
  if(eigen_values[0] == eigen_values[1])
  {
    // isotropic case (infinite number of directions)
    // by default: assemble a horizontal line that goes
    // through the centroid by default.
    line = Line(c,Direction(1,0));
  }
  else
  {
    Direction direction(eigen_vectors[0],eigen_vectors[1]);
    line = Line(c,direction);
  }
}


// this one deduces the kernel from the point in container
template < typename InputIterator, 
           typename Line>
void 
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               Line& line)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel                    K;
  return CGAL::linear_least_squares_fitting_2(begin,end,line,K());
}

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_3


