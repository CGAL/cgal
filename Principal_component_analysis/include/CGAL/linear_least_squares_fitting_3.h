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
// Author(s)     : Pierre Alliez and Sylvain Pion

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_3_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_3_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/centroid.h>
#include <CGAL/eigen.h>

#include <iterator>
#include <list>
#include <string>

CGAL_BEGIN_NAMESPACE

namespace CGALi {


// assemble covariance matrix from a point set 
template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond, 
                             typename K::FT covariance[6], // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& k,                   // kernel
                             const typename K::Point_3*)   // used for indirection
{
  typedef typename K::FT       FT;
  typedef typename K::Point_3  Point;
  typedef typename K::Vector_3 Vector;

  // Matrix numbering:
  // 0          
  // 1 2
  // 3 4 5          
  covariance[0] = covariance[1] = covariance[2] = 
  covariance[3] = covariance[4] = covariance[5] = (FT)0.0;
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Point& p = *it;
    Vector d = p - c;
    covariance[0] += d.x() * d.x();
    covariance[1] += d.x() * d.y();
    covariance[2] += d.y() * d.y();
    covariance[3] += d.x() * d.z();
    covariance[4] += d.y() * d.z();
    covariance[5] += d.z() * d.z();
  }
}

// assemble covariance matrix from a triangle set 
template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond, 
                             typename K::FT covariance[6], // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& k,                   // kernel
                             const typename K::Triangle_3*)// used for indirection
{
  typedef typename K::FT       FT;
  typedef typename K::Point_3  Point;
  typedef typename K::Vector_3 Vector;
  typedef typename K::Triangle_3 Triangle;

  // Matrix numbering:
  // 0          
  // 1 2
  // 3 4 5          
  covariance[0] = covariance[1] = covariance[2] = 
  covariance[3] = covariance[4] = covariance[5] = (FT)0.0;
  FT sum_areas = 0.0;
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Triangle& triangle = *it;
    FT area = std::sqrt(triangle.squared_area());
    Point c_t = K().construct_centroid_3_object()(triangle[0],triangle[1],triangle[2]);
    sum_areas += area;

    // e1 = ab, e2 = ac
    Vector e1 = Vector(triangle[0],triangle[1]);
    Vector e2 = Vector(triangle[0],triangle[2]);

    FT c1 = (FT)(2.0 * area * 10.0 / 72.0);
    FT c2 = (FT)(2.0 * area *  7.0 / 72.0);
        
    covariance[0] += c1*(e1[0]*e1[0] + e2[0]*e2[0]) + (FT)2.0*c2*e1[0]*e2[0];
    covariance[1] += c1*(e1[1]*e1[0] + e2[1]*e2[0]) + c2*(e1[1]*e2[0] + e1[0]*e2[1]);
    covariance[2] += c1*(e1[1]*e1[1] + e2[1]*e2[1]) + (FT)2.0*c2*e1[1]*e2[1];
    covariance[3] += c1*(e1[2]*e1[0] + e2[2]*e2[0]) + c2*(e1[2]*e2[0] + e1[0]*e2[2]);
    covariance[4] += c1*(e1[2]*e1[1] + e2[2]*e2[1]) + c2*(e1[2]*e2[1] + e1[1]*e2[2]);
    covariance[5] += c1*(e1[2]*e1[2] + e2[2]*e2[2]) + (FT)2.0*c2*e1[2]*e2[2];
    
    // add area(t) c(t) * transpose(c(t))
    covariance[0] += area * c_t.x() * c_t.x();
    covariance[1] += area * c_t.y() * c_t.x();
    covariance[2] += area * c_t.y() * c_t.y();
    covariance[3] += area * c_t.z() * c_t.x();
    covariance[4] += area * c_t.z() * c_t.y();
    covariance[5] += area * c_t.z() * c_t.z();
  }

  // remove sum(area) * (c * transpose(c))
  covariance[0] -= sum_areas * c.x() * c.x();
  covariance[1] -= sum_areas * c.y() * c.x();
  covariance[2] -= sum_areas * c.y() * c.y();
  covariance[3] -= sum_areas * c.z() * c.x();
  covariance[4] -= sum_areas * c.z() * c.y();
  covariance[5] -= sum_areas * c.z() * c.z();
}


// compute the eigen values and vectors of the covariance 
// matrix and deduces the best linear fitting plane.
// returns fitting quality
template < typename K >
typename K::FT
fitting_plane_3(const typename K::FT covariance[6], // covariance matrix
                const typename K::Point_3& c,       // centroid
                typename K::Plane_3& plane,         // best fit plane
                const K& k)                         // kernel
{
  typedef typename K::FT       FT;
  typedef typename K::Point_3  Point;
  typedef typename K::Plane_3  Plane;
  typedef typename K::Vector_3 Vector;

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  FT eigen_values[3];
  FT eigen_vectors[9];
  eigen_symmetric<FT>(covariance,3,eigen_vectors,eigen_values);

  // check unicity and build fitting line accordingly
  if(eigen_values[0] != eigen_values[1] &&
     eigen_values[0] != eigen_values[2])
  {
    // regular case
    Vector normal(eigen_vectors[6],
                  eigen_vectors[7],
                  eigen_vectors[8]);
    plane = Plane(c,normal);
    return (FT)1.0 - eigen_values[2] / eigen_values[0];
  } // end regular case
  else
  {
    // isotropic case (infinite number of directions)
    // by default: assemble a horizontal plane that goes
    // through the centroid.
    plane = Plane(c,Vector(0.0,0.0,1.0));
    return (FT)0.0;
  } 
}

// compute the eigen values and vectors of the covariance 
// matrix and deduces the best linear fitting line
// (this is an internal function)
// returns fitting quality
template < typename K >
typename K::FT
fitting_line_3(const typename K::FT covariance[6], // covariance matrix
               const typename K::Point_3& c,       // centroid
               typename K::Line_3& line,           // best fit line
               const K& k)                         // kernel
{
  typedef typename K::FT       FT;
  typedef typename K::Point_3  Point;
  typedef typename K::Line_3   Line;
  typedef typename K::Vector_3 Vector;

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in descending order, 
  // eigen vectors are sorted in accordance.
  FT eigen_values[3];
  FT eigen_vectors[9];
  eigen_symmetric<FT>(covariance,3,eigen_vectors,eigen_values);

  // check unicity and build fitting line accordingly
  if(eigen_values[0] != eigen_values[1])
  {
    // regular case
    Vector direction(eigen_vectors[0],eigen_vectors[1],eigen_vectors[2]);
    line = Line(c,direction);
    return (FT)1.0 - eigen_values[1] / eigen_values[0];
  } // end regular case
  else
  {
    // isotropic case (infinite number of directions)
    // by default: assemble a horizontal plane that goes
    // through the centroid.
    line = Line(c,Vector(1.0,0.0,0.0));
    return (FT)0.0;
  } 
}



// fits a plane to a 3D point set
// returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best (zero variance orthogonally to the fitting line)
//  0 is worst (isotropic case, returns a plane with default direction)
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane, // best fit plane
                               typename K::Point_3& c,     // centroid
                               const K& k,                 // kernel
                               const typename K::Point_3*) // used for indirection
{
  typedef typename K::FT       FT;
  typedef typename K::Point_3  Point;
  typedef typename K::Plane_3  Plane;
  typedef typename K::Vector_3 Vector;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K());

  // assemble covariance matrix
  FT covariance[6];
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Point*) NULL);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k);
} // end fit plane to point set

// fits a line to a 3D point set
// returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best (zero variance orthogonally to the fitting line)
//  0 is worst (isotropic case, returns a line along x axis)
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,  // best fit line
                               typename K::Point_3& c,    // centroid
                               const K& k,                // kernel
                               const typename K::Point_3*)// used for indirection
{
  typedef typename K::FT       FT;
  typedef typename K::Point_3  Point;
  typedef typename K::Line_3   Line;
  typedef typename K::Vector_3 Vector;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K());

  // assemble covariance matrix
  FT covariance[6];
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Point*) NULL);

  // compute fitting line
  return fitting_line_3(covariance,c,line,k);
} // end fit line to point set

// fits a plane to a 3D triangle set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,   // best fit plane
                               typename K::Point_3& c,       // centroid
                               const K& k,                   // kernel
                               const typename K::Triangle_3*)// used for indirection
{
  typedef typename K::FT          FT;
  typedef typename K::Point_3     Point;
  typedef typename K::Vector_3    Vector;
  typedef typename K::Triangle_3  Triangle;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K());

  // assemble covariance matrix
  FT covariance[6];
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Triangle*) NULL);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k);
} // end linear_least_squares_fitting_3

// fits a line to a 3D triangle set
template < typename InputIterator, 
           typename K >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,     // best fit line
                               typename K::Point_3& c,       // centroid
                               const K& k,                   // kernel
                               const typename K::Triangle_3*)// used for indirection
{
  typedef typename K::FT          FT;
  typedef typename K::Point_3     Point;
  typedef typename K::Vector_3    Vector;
  typedef typename K::Triangle_3  Triangle;
  typedef typename K::Line_3      Line;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K());

  // assemble covariance matrix
  FT covariance[6];
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Triangle*) NULL);

  // compute fitting plane
  return fitting_line_3(covariance,c,line,k);
} // end linear_least_squares_fitting_3

} // end namespace CGALi

template < typename InputIterator, 
           typename K >
inline
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,
                               typename K::Point_3& centroid,
                               const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return CGALi::linear_least_squares_fitting_3(first, beyond, plane,
                                               centroid, k, (Value_type*) NULL);
}

template < typename InputIterator, 
           typename K >
inline
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,
                               typename K::Point_3& centroid,
                               const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  return CGALi::linear_least_squares_fitting_3(first, beyond, line,
                                               centroid, k, (Value_type*) NULL);
}

template < typename InputIterator, 
           typename K >
inline
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Plane_3& plane,
                               const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typename K::Point_3 centroid;
  return CGALi::linear_least_squares_fitting_3(first, beyond, plane,
                                               centroid, k,(Value_type*) NULL);
}

template < typename InputIterator, 
           typename K >
inline
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               typename K::Line_3& line,
                               const K& k)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typename K::Point_3 centroid;
  return CGALi::linear_least_squares_fitting_3(first, beyond, line,
                                               centroid, k,(Value_type*) NULL);
}

// deduces the kernel from the points in container.
template < typename InputIterator, 
           typename Object,
           typename Point>
inline
typename Kernel_traits<Object>::Kernel::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond, 
                               Object& object,
                               Point& centroid)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel K;
  return CGAL::linear_least_squares_fitting_3(first,beyond,object,centroid,K());
}

// does not return the centroid and deduces the kernel as well.
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

CGAL_END_NAMESPACE

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_3_H
