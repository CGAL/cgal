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
#include <CGAL/eigen_2.h>
#include <CGAL/centroid.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// fit a line to a 2D point set
// return a fitting quality (1 - lambda_min/lambda_max):
//  1 is best (zero variance orthogonally to the fitting line)
//  0 is worst (isotropic case, return a line with default direction)
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               typename K::Line_2& line,
                               const K& k, 
			       const typename K::Point_2*)
{
  typedef typename K::FT          FT;
  typedef typename K::Point_2     Point;
  typedef typename K::Vector_2    Vector;
  typedef typename K::Line_2      Line;

  // precondition: at least one element in the container.
  CGAL_precondition(begin != end);

  // compute centroid
  Point c = centroid(begin,end,K());

  // assemble covariance matrix as a
  // semi-definite matrix. 
  // Matrix numbering:
  // 0
  // 1 2
  FT covariance[3] = {0,0,0};
  std::pair<FT,FT> eigen_values;
  std::pair<Vector,Vector> eigen_vectors;
  for(InputIterator it = begin;
		    it != end;
		    it++)
  {
    const Point& p = *it;
    Vector d = p - c;
    covariance[0] += d.x() * d.x();
    covariance[1] += d.x() * d.y();
    covariance[2] += d.y() * d.y();
  }

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  // TODO: use explicit formula instead.
  eigen_symmetric_2<K>(covariance,eigen_vectors,eigen_values);

  // assert eigen values are positives
  CGAL_assertion(eigen_values.first >= 0 && 
                       eigen_values.second >= 0);

  // check unicity and build fitting line accordingly
  if(eigen_values.first != eigen_values.second)
  {
    // regular case
    line = Line(c,eigen_vectors.first);
    return (FT)1.0 - eigen_values.second / eigen_values.first;
  } // end regular case
  else
  {
    // isotropic case (infinite number of directions)
    // by default: assemble a line that goes through 
    // the centroid and with a default vector.
    line = Line(c,Vector(1,0));
    return (FT)0.0;
  } // end isotropic case
} // end linear_least_squares_fitting_2

// fit a line to a 2D triangle set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               typename K::Line_2& line,
                               const K& k, const typename K::Triangle_2*)
{
  typedef typename K::FT          FT;
  typedef typename K::Point_2     Point;
  typedef typename K::Vector_2    Vector;
  typedef typename K::Triangle_2  Triangle;
  typedef typename K::Line_2      Line;

  // precondition: at least one element in the container.
  CGAL_precondition(begin != end);

  // compute centroid
  Point c = centroid(begin,end,K());

  // assemble covariance matrix as a
  // semi-definite matrix. 
  // Matrix numbering:
  // 0          
  // 1 2
  FT covariance[3] = {0,0,0};
  FT eigen_values[2] = {0,0};
  // these should be typed Vector_2
  // but it requires working more on eigen.h
  FT eigen_vectors[4] = {0,0,0,0};
  FT sum_areas = 0;
  for(InputIterator it = begin;
		    it != end;
		    it++)
  {
    const Triangle& triangle = *it;
    FT area = std::abs(triangle.area());
    Point g = centroid(triangle); // local centroid
    sum_areas += area;

    // e1 = ab, e2 = ac
    Vector e1 = triangle[1] - triangle[0];
    Vector e2 = triangle[2] - triangle[0];

    FT coef1 = 2.0 * area * 10.0/72.0;
    FT coef2 = 2.0 * area * 7.0/72.0;
        
    covariance[0] += coef1*(e1[0]*e1[0] + e2[0]*e2[0]) + 2.0*coef2*e1[0]*e2[0];
    covariance[1] += coef1*(e1[1]*e1[0] + e2[1]*e2[0]) + coef2*(e1[1]*e2[0] + e1[0]*e2[1]);
    covariance[2] += coef1*(e1[1]*e1[1] + e2[1]*e2[1]) + 2.0*coef2*e1[1]*e2[1];
    
    // add area(t) g(t)*transpose(g(t))
    covariance[0] += area * g.x() * g.x();
    covariance[1] += area * g.y() * g.x();
    covariance[2] += area * g.y() * g.y();
  }

  // remove sum_t(area) * (c * transpose(c))
  covariance[0] -= sum_areas * c.x() * c.x();
  covariance[1] -= sum_areas * c.y() * c.x();
  covariance[2] -= sum_areas * c.y() * c.y();
  
  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  // TODO: use explicit formula instead.
  eigen_semi_definite_symmetric(covariance,2,eigen_vectors,eigen_values);

  // assert eigen values are positives
  CGAL_assertion(eigen_values[0] >= 0 && 
                       eigen_values[1] >= 0);

  // check unicity and build fitting line accordingly
  if(eigen_values[0] != eigen_values[1])
  {
    // regular case
    line = Line(c,Vector(eigen_vectors[0],eigen_vectors[1]));
    return (FT)1.0 - eigen_values[1] / eigen_values[0];
  } // end regular case
  else
  {
    // isotropic case (infinite number of directions)
    // by default: assemble a line that goes through 
    // the centroid and with a default horizontal direction.
    line = Line(c,Vector(1,0));
    return (FT)0.0;
  } // end isotropic case
} // end linear_least_squares_fitting_2

} // namespace CGALi


template < typename InputIterator, 
           typename K >
inline
typename K::FT
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               typename K::Line_2& line,
                               const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type
    Value_type;
  return CGALi::linear_least_squares_fitting_2(begin, end, line, k,
					       (Value_type*) NULL);
}


// this one deduces the kernel from the point in container
template < typename InputIterator, 
           typename Line >
inline
typename Kernel_traits<Line>::Kernel::FT
linear_least_squares_fitting_2(InputIterator begin,
                               InputIterator end, 
                               Line& line)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel               K;
  return CGAL::linear_least_squares_fitting_2(begin,end,line,K());
}

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_2


