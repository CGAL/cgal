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
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_UTIL_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_UTIL_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/Linear_algebraCd.h>
#include <CGAL/Dimension.h>

namespace CGAL {

namespace internal {

// Initialize a matrix in n dimension by an array or numbers
template <typename K>
typename CGAL::Linear_algebraCd<typename K::FT>::Matrix
init_matrix(const int n,
            typename K::FT entries[])
{
  CGAL_assertion(n > 1); // dimension > 1
  typedef typename CGAL::Linear_algebraCd<typename K::FT>::Matrix Matrix;

  Matrix m(n);
  int i,j;
  for(i = 0; i < n; i++) 
    for(j = 0; j < n; j++)
      m[i][j] = entries[i*n+j];

  return m;
} // end initialization of matrix

// assemble covariance matrix from a point set 
template < typename InputIterator,
           typename K,
	   typename DiagonalizeTraits >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond, 
                             typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& k,                    // kernel
                             const typename K::Point_3*,   // used for indirection
                             const CGAL::Dimension_tag<0>&,
			     const DiagonalizeTraits&)
{
  typedef typename K::FT       FT;
  typedef typename K::Point_3  Point;
  typedef typename K::Vector_3 Vector;

  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5          
  covariance[0] = covariance[1] = covariance[2] = 
  covariance[3] = covariance[4] = covariance[5] = (FT)0.0;
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Point& p = *it;
    Vector d = k.construct_vector_3_object()(c,p);
    covariance[0] += d.x() * d.x();
    covariance[1] += d.x() * d.y();
    covariance[2] += d.x() * d.z();
    covariance[3] += d.y() * d.y();
    covariance[4] += d.y() * d.z();
    covariance[5] += d.z() * d.z();
  }
}

// assemble covariance matrix from a triangle set 
template < typename InputIterator,
           typename K,
	   typename DiagonalizeTraits >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
			     typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K&,                    // kernel
                             const typename K::Triangle_3*,// used for indirection
                             const CGAL::Dimension_tag<2>&,
			     const DiagonalizeTraits&)
{
  typedef typename K::FT          FT;
  typedef typename K::Triangle_3  Triangle;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // assemble covariance matrix as a semi-definite matrix. 
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5          
  //Final combined covariance matrix for all triangles and their combined mass
  FT mass = 0.0;

  // assemble 2nd order moment about the origin.  
  FT temp[9] = {1.0/12.0, 1.0/24.0, 1.0/24.0,
                1.0/24.0, 1.0/12.0, 1.0/24.0,
                1.0/24.0, 1.0/24.0, 1.0/12.0};
  Matrix moment = init_matrix<K>(3,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each triangle, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Triangle& t = *it;

    // defined for convenience.
    FT delta[9] = {t[0].x(), t[1].x(), t[2].x(), 
                   t[0].y(), t[1].y(), t[2].y(),
                   t[0].z(), t[1].z(), t[2].z()};
    Matrix transformation = init_matrix<K>(3,delta);
    FT area = std::sqrt(t.squared_area());

		// skip zero measure primitives
    if(area == (FT)0.0)
			continue;

    // Find the 2nd order moment for the triangle wrt to the origin by an affine transformation.
    
    // Transform the standard 2nd order moment using the transformation matrix
    transformation = 2 * area * transformation * moment * LA::transpose(transformation);
    
    // and add to covariance matrix
    covariance[0] += transformation[0][0];
    covariance[1] += transformation[1][0];
    covariance[2] += transformation[2][0];
    covariance[3] += transformation[1][1];
    covariance[4] += transformation[2][1];
    covariance[5] += transformation[2][2];

    mass += area;
  }

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += mass * (-1.0 * c.x() * c.x());
  covariance[1] += mass * (-1.0 * c.x() * c.y());
  covariance[2] += mass * (-1.0 * c.z() * c.x());
  covariance[3] += mass * (-1.0 * c.y() * c.y());
  covariance[4] += mass * (-1.0 * c.z() * c.y());
  covariance[5] += mass * (-1.0 * c.z() * c.z());

}

// assemble covariance matrix from a cuboid set 
template < typename InputIterator,
           typename K,
	   typename DiagonalizeTraits >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
			     typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& ,                    // kernel
                             const typename K::Iso_cuboid_3*,// used for indirection
                             const CGAL::Dimension_tag<3>&,
			     const DiagonalizeTraits&)
{
  typedef typename K::FT          FT;
  typedef typename K::Iso_cuboid_3    Iso_cuboid;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // assemble covariance matrix as a semi-definite matrix. 
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5          
  // final combined covariance matrix for all cuboids and their combined mass
  FT mass = (FT)0.0;

  // assemble 2nd order moment about the origin.  
  FT temp[9] = {(FT)(1.0/3.0), (FT)(1.0/4.0), (FT)(1.0/4.0),
								(FT)(1.0/4.0), (FT)(1.0/3.0), (FT)(1.0/4.0),
                (FT)(1.0/4.0), (FT)(1.0/4.0), (FT)(1.0/3.0)};
  Matrix moment = init_matrix<K>(3,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each cuboid, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Iso_cuboid& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT x0 = t[0].x();
    FT y0 = t[0].y();
    FT z0 = t[0].z();
    FT delta[9] = {t[1].x()-x0, t[3].x()-x0, t[5].x()-x0, 
                   t[1].y()-y0, t[3].y()-y0, t[5].y()-y0,
                   t[1].z()-z0, t[3].z()-z0, t[5].z()-z0};
    Matrix transformation = init_matrix<K>(3,delta);
    FT volume = t.volume();

		// skip zero measure primitives
    if(volume == (FT)0.0)
			continue;

    // Find the 2nd order moment for the cuboid wrt to the origin by an affine transformation.
    
    // Transform the standard 2nd order moment using the transformation matrix
    transformation = volume * transformation * moment * LA::transpose(transformation);
    
    // Translate the 2nd order moment to the minimum corner (x0,y0,z0) of the cuboid.
    FT xav0 = (delta[0] + delta[1] + delta[2])/4.0;
    FT yav0 = (delta[3] + delta[4] + delta[5])/4.0;
    FT zav0 = (delta[6] + delta[7] + delta[8])/4.0;

    // and add to covariance matrix
    covariance[0] += transformation[0][0] + volume * (2*x0*xav0 + x0*x0);
    covariance[1] += transformation[1][0] + volume * (xav0*y0 + yav0*x0 + x0*y0);
    covariance[2] += transformation[2][0] + volume * (x0*zav0 + xav0*z0 + x0*z0);
    covariance[3] += transformation[1][1] + volume * (2*y0*yav0 + y0*y0);
    covariance[4] += transformation[2][1] + volume * (yav0*z0 + y0*zav0 + z0*y0);
    covariance[5] += transformation[2][2] + volume * (2*zav0*z0 + z0*z0);

    mass += volume;
  }

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += mass * (- c.x() * c.x());
  covariance[1] += mass * (- c.x() * c.y());
  covariance[2] += mass * (- c.z() * c.x());
  covariance[3] += mass * (- c.y() * c.y());
  covariance[4] += mass * (- c.z() * c.y());
  covariance[5] += mass * (- c.z() * c.z());
}

// assemble covariance matrix from a cuboid set 
template < typename InputIterator,
           typename K,
	   typename DiagonalizeTraits >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond, 
			     typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& ,                    // kernel
                             const typename K::Iso_cuboid_3*,// used for indirection
                             const CGAL::Dimension_tag<2>&,
			     const DiagonalizeTraits&)
{
  typedef typename K::FT FT;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // assemble covariance matrix as a semi-definite matrix. 
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5          
  //Final combined covariance matrix for all cuboids and their combined mass
  FT mass = (FT)0.0;

  // assemble 2nd order moment about the origin.  
  FT temp[9] = {(FT)(7.0/3.0), (FT)1.5,       (FT)1.5,
                (FT)1.5,       (FT)(7.0/3.0), (FT)1.5,
                (FT)1.5,       (FT)1.5,       (FT)(7.0/3.0)};
  Matrix moment = init_matrix<K>(3,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each cuboid, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Iso_cuboid& t = *it;

    // defined for convenience.
    FT x0 = t[0].x();
    FT y0 = t[0].y();
    FT z0 = t[0].z();
    FT delta[9] = {t[1].x()-x0, t[3].x()-x0, t[5].x()-x0, 
                   t[1].y()-y0, t[3].y()-y0, t[5].y()-y0,
                   t[1].z()-z0, t[3].z()-z0, t[5].z()-z0};
    Matrix transformation = init_matrix<K>(3,delta);
    FT area = std::pow(delta[0]*delta[0] + delta[3]*delta[3] +
                  delta[6]*delta[6],1/3.0)*std::pow(delta[1]*delta[1] +
                  delta[4]*delta[4] + delta[7]*delta[7],1/3.0)*2 +
                  std::pow(delta[0]*delta[0] + delta[3]*delta[3] +
                  delta[6]*delta[6],1/3.0)*std::pow(delta[2]*delta[2] +
                  delta[5]*delta[5] + delta[8]*delta[8],1/3.0)*2 +
                  std::pow(delta[1]*delta[1] + delta[4]*delta[4] +
                  delta[7]*delta[7],1/3.0)*std::pow(delta[2]*delta[2] +
                  delta[5]*delta[5] + delta[8]*delta[8],1/3.0)*2;

		// skip zero measure primitives
    if(area == (FT)0.0)
			continue;

    // Find the 2nd order moment for the cuboid wrt to the origin by an affine transformation.
    
    // Transform the standard 2nd order moment using the transformation matrix
    transformation = area * transformation * moment * LA::transpose(transformation);
    
    // Translate the 2nd order moment to the minimum corner (x0,y0,z0) of the cuboid.
    FT xav0 = (delta[0] + delta[1] + delta[2])/4.0;
    FT yav0 = (delta[3] + delta[4] + delta[5])/4.0;
    FT zav0 = (delta[6] + delta[7] + delta[8])/4.0;

    // and add to covariance matrix
    covariance[0] += transformation[0][0] + area * (2*x0*xav0 + x0*x0);
    covariance[1] += transformation[1][0] + area * (xav0*y0 + yav0*x0 + x0*y0);
    covariance[2] += transformation[2][0] + area * (x0*zav0 + xav0*z0 + x0*z0);
    covariance[3] += transformation[1][1] + area * (2*y0*yav0 + y0*y0);
    covariance[4] += transformation[2][1] + area * (yav0*z0 + y0*zav0 + z0*y0);
    covariance[5] += transformation[2][2] + area * (2*zav0*z0 + z0*z0);

    mass += area;
  }

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += mass * (-1.0 * c.x() * c.x());
  covariance[1] += mass * (-1.0 * c.x() * c.y());
  covariance[2] += mass * (-1.0 * c.z() * c.x());
  covariance[3] += mass * (-1.0 * c.y() * c.y());
  covariance[4] += mass * (-1.0 * c.z() * c.y());
  covariance[5] += mass * (-1.0 * c.z() * c.z());

}

// assemble covariance matrix from a sphere set 
template < typename InputIterator,
           typename K,
	   typename DiagonalizeTraits >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond, 
			     typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K&,                     // kernel
                             const typename K::Sphere_3*,  // used for indirection
                             const CGAL::Dimension_tag<3>&,
			     const DiagonalizeTraits&)
{
  typedef typename K::FT          FT;
  typedef typename K::Sphere_3  Sphere;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // assemble covariance matrix as a semi-definite matrix. 
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5          
  //Final combined covariance matrix for all spheres and their combined mass
  FT mass = 0.0;

  // assemble 2nd order moment about the origin.  
  FT temp[9] = {4.0/15.0, 0.0,      0.0,
                0.0,      4.0/15.0, 0.0,
                0.0,      0.0,      4.0/15.0};
  Matrix moment = init_matrix<K>(3,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each sphere, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Sphere& t = *it;

    // defined for convenience.
    FT radius = std::sqrt(t.squared_radius());
    FT delta[9] = {radius, 0.0, 0.0, 
                   0.0, radius, 0.0,
                   0.0, 0.0, radius};
    Matrix transformation = init_matrix<K>(3,delta);
    FT volume = (FT)(4.0/3.0) * radius * t.squared_radius();

		// skip zero measure primitives
    if(volume == (FT)0.0)
			continue;

    // Find the 2nd order moment for the sphere wrt to the origin by an affine transformation.
    
    // Transform the standard 2nd order moment using the transformation matrix
    transformation = (3.0/4.0) * volume * transformation * moment * LA::transpose(transformation);
    
    // Translate the 2nd order moment to the center of the sphere.
    FT x0 = t.center().x();
    FT y0 = t.center().y();
    FT z0 = t.center().z();

    // and add to covariance matrix
    covariance[0] += transformation[0][0] + volume * x0*x0;
    covariance[1] += transformation[1][0] + volume * x0*y0;
    covariance[2] += transformation[2][0] + volume * x0*z0;
    covariance[3] += transformation[1][1] + volume * y0*y0;
    covariance[4] += transformation[2][1] + volume * z0*y0;
    covariance[5] += transformation[2][2] + volume * z0*z0;

    mass += volume;
  }

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += mass * (-1.0 * c.x() * c.x());
  covariance[1] += mass * (-1.0 * c.x() * c.y());
  covariance[2] += mass * (-1.0 * c.z() * c.x());
  covariance[3] += mass * (-1.0 * c.y() * c.y());
  covariance[4] += mass * (-1.0 * c.z() * c.y());
  covariance[5] += mass * (-1.0 * c.z() * c.z());

}
// assemble covariance matrix from a sphere set 
template < typename InputIterator,
           typename K,
	   typename DiagonalizeTraits >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond, 
			     typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K&,                     // kernel
                             const typename K::Sphere_3*,  // used for indirection
                             const CGAL::Dimension_tag<2>&,
			     const DiagonalizeTraits&)
{
  typedef typename K::FT          FT;
  typedef typename K::Sphere_3  Sphere;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // assemble covariance matrix as a semi-definite matrix. 
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5          
  //Final combined covariance matrix for all spheres and their combined mass
  FT mass = 0.0;

  // assemble 2nd order moment about the origin.  
  FT temp[9] = {4.0/3.0, 0.0,     0.0,
                0.0,     4.0/3.0, 0.0,
                0.0,     0.0,     4.0/3.0};
  Matrix moment = init_matrix<K>(3,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each sphere, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Sphere& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT radius = std::sqrt(t.squared_radius());
    FT delta[9] = {radius, 0.0,    0.0, 
                   0.0,    radius, 0.0,
                   0.0,    0.0,    radius};
    Matrix transformation = init_matrix<K>(3,delta);
    FT area = (FT)4.0 * t.squared_radius();

		// skip zero measure primitives
    if(area == (FT)0.0)
			continue;

    // Find the 2nd order moment for the sphere wrt to the origin by an affine transformation.
    
    // Transform the standard 2nd order moment using the transformation matrix
    transformation = (1.0/4.0) * area * transformation * moment * LA::transpose(transformation);
    
    // Translate the 2nd order moment to the center of the sphere.
    FT x0 = t.center().x();
    FT y0 = t.center().y();
    FT z0 = t.center().z();

    // and add to covariance matrix
    covariance[0] += transformation[0][0] + area * x0*x0;
    covariance[1] += transformation[1][0] + area * x0*y0;
    covariance[2] += transformation[2][0] + area * x0*z0;
    covariance[3] += transformation[1][1] + area * y0*y0;
    covariance[4] += transformation[2][1] + area * z0*y0;
    covariance[5] += transformation[2][2] + area * z0*z0;

    mass += area;
  }

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += mass * (-1.0 * c.x() * c.x());
  covariance[1] += mass * (-1.0 * c.x() * c.y());
  covariance[2] += mass * (-1.0 * c.z() * c.x());
  covariance[3] += mass * (-1.0 * c.y() * c.y());
  covariance[4] += mass * (-1.0 * c.z() * c.y());
  covariance[5] += mass * (-1.0 * c.z() * c.z());

}

// assemble covariance matrix from a tetrahedron set 
template < typename InputIterator,
           typename K,
	   typename DiagonalizeTraits >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond, 
			     typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& ,                    // kernel
                             const typename K::Tetrahedron_3*,// used for indirection
                             const CGAL::Dimension_tag<3>&,
			     const DiagonalizeTraits&)
{
  typedef typename K::FT          FT;
  typedef typename K::Tetrahedron_3  Tetrahedron;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // assemble covariance matrix as a semi-definite matrix. 
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5          
  //Final combined covariance matrix for all tetrahedrons and their combined mass
  FT mass = 0.0;

  // assemble 2nd order moment about the origin.  
  FT temp[9] = {1.0/60.0,  1.0/120.0, 1.0/120.0,
                1.0/120.0, 1.0/60.0,  1.0/120.0,
                1.0/120.0, 1.0/120.0, 1.0/60.0};
  Matrix moment = init_matrix<K>(3,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each tetrahedron, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Tetrahedron& t = *it;

    // defined for convenience.
    FT x0 = t[0].x();
    FT y0 = t[0].y();
    FT z0 = t[0].z();

    FT delta[9] = {t[1].x()-x0, t[2].x()-x0, t[3].x()-x0, 
                   t[1].y()-y0, t[2].y()-y0, t[3].y()-y0,
                   t[1].z()-z0, t[2].z()-z0, t[3].z()-z0};
    Matrix transformation = init_matrix<K>(3,delta);
    FT volume = t.volume();

		// skip zero measure primitives
    if(volume == (FT)0.0)
			continue;

    // Find the 2nd order moment for the tetrahedron wrt to the origin by an affine transformation.
    
    // Transform the standard 2nd order moment using the transformation matrix
    transformation = 6 * volume * transformation * moment * LA::transpose(transformation);
    
    // Translate the 2nd order moment to the center of the tetrahedron.
    FT xav0 = (delta[0]+delta[1]+delta[2])/4.0;
    FT yav0 = (delta[3]+delta[4]+delta[5])/4.0;
    FT zav0 = (delta[6]+delta[7]+delta[8])/4.0;

    // and add to covariance matrix
    covariance[0] += transformation[0][0] + volume * (2*x0*xav0 + x0*x0);
    covariance[1] += transformation[1][0] + volume * (xav0*y0 + yav0*x0 + x0*y0);
    covariance[2] += transformation[2][0] + volume * (x0*zav0 + xav0*z0 + x0*z0);
    covariance[3] += transformation[1][1] + volume * (2*y0*yav0 + y0*y0);
    covariance[4] += transformation[2][1] + volume * (yav0*z0 + y0*zav0 + z0*y0);
    covariance[5] += transformation[2][2] + volume * (2*zav0*z0 + z0*z0);

    mass += volume;
  }

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += mass * (-1.0 * c.x() * c.x());
  covariance[1] += mass * (-1.0 * c.x() * c.y());
  covariance[2] += mass * (-1.0 * c.z() * c.x());
  covariance[3] += mass * (-1.0 * c.y() * c.y());
  covariance[4] += mass * (-1.0 * c.z() * c.y());
  covariance[5] += mass * (-1.0 * c.z() * c.z());
}

// assemble covariance matrix from a segment set 
template < typename InputIterator,
           typename K,
	   typename DiagonalizeTraits >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond, 
			     typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& ,                    // kernel
                             const typename K::Segment_3*,// used for indirection
                             const CGAL::Dimension_tag<1>&,
			     const DiagonalizeTraits&)
{
  typedef typename K::FT          FT;
  typedef typename K::Segment_3  Segment;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // assemble covariance matrix as a semi-definite matrix. 
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5          
  //Final combined covariance matrix for all segments and their combined mass
  FT mass = 0.0;

  // assemble 2nd order moment about the origin.  
  FT temp[9] = {1.0, 0.5, 0.0,
                0.5, 1.0, 0.0,
                0.0, 0.0, 0.0};
  Matrix moment = (FT)(1.0/3.0) * init_matrix<K>(3,temp);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each segment, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Segment& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT delta[9] = {t[0].x(), t[1].x(), 0.0, 
       t[0].y(), t[1].y(), 0.0,
                   t[0].z(), t[1].z(), 1.0};
    Matrix transformation = init_matrix<K>(3,delta);
    FT length = std::sqrt(t.squared_length());

		// skip zero measure primitives
    if(length == (FT)0.0)
			continue;

    // Find the 2nd order moment for the segment wrt to the origin by an affine transformation.
    
    // Transform the standard 2nd order moment using the transformation matrix
    transformation = length * transformation * moment * LA::transpose(transformation);

    // and add to covariance matrix
    covariance[0] += transformation[0][0];
    covariance[1] += transformation[1][0];
    covariance[2] += transformation[2][0];
    covariance[3] += transformation[1][1];
    covariance[4] += transformation[2][1];
    covariance[5] += transformation[2][2];

    mass += length;
  }

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += mass * (-1.0 * c.x() * c.x());
  covariance[1] += mass * (-1.0 * c.x() * c.y());
  covariance[2] += mass * (-1.0 * c.z() * c.x());
  covariance[3] += mass * (-1.0 * c.y() * c.y());
  covariance[4] += mass * (-1.0 * c.z() * c.y());
  covariance[5] += mass * (-1.0 * c.z() * c.z());

}


// compute the eigen values and vectors of the covariance 
// matrix and deduces the best linear fitting plane.
// returns fitting quality
template < typename K, typename DiagonalizeTraits >
typename K::FT
fitting_plane_3(typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
                const typename K::Point_3& c,       // centroid
                typename K::Plane_3& plane,         // best fit plane
                const K&,                           // kernel
		const DiagonalizeTraits& )                 // Diagonalize traits
{
  typedef typename K::FT       FT;
  typedef typename K::Plane_3  Plane;
  typedef typename K::Vector_3 Vector;

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in ascending order, 
  // eigen vectors are sorted in accordance.
  typename DiagonalizeTraits::Vector eigen_values = {{ 0. , 0., 0. }};
  typename DiagonalizeTraits::Matrix eigen_vectors = {{ 0., 0., 0.,
							0., 0., 0.,
							0., 0., 0. }};
  DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
    (covariance, eigen_values, eigen_vectors);

  // degenerate case 
  if(eigen_values[0] == eigen_values[1] && 
     eigen_values[1] == eigen_values[2])
  {
    // assemble a default horizontal plane that goes
    // through the centroid.
    plane = Plane(c,Vector(FT(0),FT(0),FT(1)));
    return FT(0);
  } 
  else // regular and line case
  {
    Vector normal(eigen_vectors[0],
                  eigen_vectors[1],
                  eigen_vectors[2]);
    plane = Plane(c,normal);
    return FT(1) - eigen_values[0] / eigen_values[1];
  } // end regular case
}

// compute the eigen values and vectors of the covariance 
// matrix and deduces the best linear fitting line
// (this is an internal function)
// returns fitting quality
template < typename K, typename DiagonalizeTraits >
typename K::FT
fitting_line_3(typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
               const typename K::Point_3& c,       // centroid
               typename K::Line_3& line,           // best fit line
	       const K&,                           // kernel
	       const DiagonalizeTraits& )                 // Diagonalize traits
{
  typedef typename K::FT       FT;
  typedef typename K::Line_3   Line;
  typedef typename K::Vector_3 Vector;

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in ascending order, 
  // eigen vectors are sorted in accordance.
  typename DiagonalizeTraits::Vector eigen_values = {{ 0. , 0., 0. }};
  typename DiagonalizeTraits::Matrix eigen_vectors = {{ 0., 0., 0.,
							0., 0., 0.,
							0., 0., 0. }};
  DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
    (covariance, eigen_values, eigen_vectors);

    // isotropic case (infinite number of directions)
  if(eigen_values[0] == eigen_values[1] && 
     eigen_values[0] == eigen_values[2])
  {
    // assemble a default line along x axis which goes
    // through the centroid.
    line = Line(c,Vector(FT(1),FT(0),FT(0)));
    return (FT)0.0;
  }
  else
  {
    // regular case
    Vector direction(eigen_vectors[6],eigen_vectors[7],eigen_vectors[8]);
    line = Line(c,direction);
    return (FT)1.0 - eigen_values[1] / eigen_values[2];
  } 
}

} // end namespace internal

} //namespace CGAL

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/PCA_util_Eigen.h>
#endif


#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_UTIL_H
