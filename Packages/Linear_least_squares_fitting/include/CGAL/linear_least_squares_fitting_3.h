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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_3
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_3

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <iterator>
#include <list>
#include <string>
#include <CGAL/eigen.h>

CGAL_BEGIN_NAMESPACE

// should return quality
template <class InputIterator, class Subspace, class Traits>
void linear_least_squares_fitting_3(InputIterator first,
                                    InputIterator past, 
                                    Subspace& subspace,
                                    const Traits& fitting_traits) // todo: add default traits
{
  // preconditions
  CGAL_precondition(first != past);

  linear_least_squares_fitting_3(first,past,subspace,fitting_traits);
}

// compute centroid of a point set
template <class InputIterator, class Point, class Traits>
Point centroid(InputIterator first,
               InputIterator past,
               const Traits& fitting_traits)
{
  typedef typename Traits::FT FT;
  typedef typename Traits::Vector_3 Vector;

  Vector sum = CGAL::NULL_VECTOR;
  unsigned int nb_points = 0;
  for(InputIterator it = first;
			it != past;
			it++,nb_points++)
    sum = sum + (*it-CGAL::ORIGIN);
  return (CGAL::ORIGIN + sum / (FT)nb_points);
}

template <class InputIterator, class Plane, class Traits>
void linear_least_squares_fitting_3(InputIterator first,
                                    InputIterator past, 
                                    Plane& plane,
                                    const Traits& fitting_traits)
{
  // types
  typedef typename Traits::FT FT;
  typedef typename Traits::Line_3 Line;
  typedef typename Traits::Point_3 Point;
  typedef typename Traits::Vector_3 Vector;
  typedef typename Traits::Direction_3 Direction;

  Point centroid = centroid<InputIterator,Point>(first,past,fitting_traits);

  // build covariance matrix (semi-definite)
  FT covariance[6] = {0,
                      0,0,
                      0,0,0};
  FT eigen_values[3] = {0,0,0};
  FT eigen_vectors[9] = {0,0,0,
                         0,0,0,
                         0,0,0};
  for(InputIterator it = first;
			it != past;
			it++)
  {
    const Point& p = *it;
    // numbering in covariance matrix
    // 0          
    // 1 2        
    // 3 4 5      
    FT dx = p.x()-centroid.x();
    FT dy = p.y()-centroid.y();
    FT dz = p.z()-centroid.z();
    covariance[0] += dx * dx;
    covariance[1] += dx * dy;
    covariance[2] += dy * dy;
    covariance[3] += dx * dz;
    covariance[4] += dy * dz;
    covariance[5] += dz * dz;
  }

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  eigen_semi_definite_symmetric(covariance,3,eigen_vectors,eigen_values);

  // check unicity and build fitting plane accordingly
  if(eigen_values[0] == eigen_values[1] &&
     eigen_values[1] == eigen_values[2])
  {
    // assemble a horizontal plane that goes
    // through the centroid by default.
    Direction direction(0,0,1);
    plane = Plane(centroid,direction);
  }
  else
  {
    Direction direction(eigen_vectors[6],
                        eigen_vectors[7],
                        eigen_vectors[8]);
    plane = Plane(centroid,direction);
  }
}

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_3


