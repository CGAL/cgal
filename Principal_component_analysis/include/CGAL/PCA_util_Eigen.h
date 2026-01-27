// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta and Simon Giraudot

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_UTIL_EIGEN_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_UTIL_EIGEN_H

#include <CGAL/license/Principal_component_analysis.h>

#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/Dimension.h>
#include <CGAL/Subiterator.h>

namespace CGAL {

namespace internal {

template <typename FT>
FT approximate_cbrt (const FT& x)
{
  return static_cast<FT>(std::cbrt (CGAL::to_double(x)));
}

// assemble covariance matrix from a triangle set
template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Eigen_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K&,                    // kernel
                             const typename K::Triangle_3*,// used for indirection
                             const CGAL::Dimension_tag<2>&,
                             const Eigen_diagonalize_traits<typename K::FT, 3>&)
{
  typedef typename K::FT          FT;
  typedef typename K::Triangle_3  Triangle;
  typedef typename Eigen::Matrix<FT, 3, 3> Matrix;

  // assemble covariance matrix as a semi-definite matrix.
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5
  //Final combined covariance matrix for all triangles and their combined mass
  FT mass = 0.0;

  // assemble 2nd order moment about the origin.
  Matrix moment;
  moment << FT(1.0/12.0), FT(1.0/24.0), FT(1.0/24.0),
            FT(1.0/24.0), FT(1.0/12.0), FT(1.0/24.0),
            FT(1.0/24.0), FT(1.0/24.0), FT(1.0/12.0);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each triangle, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Triangle& t = *it;

    // defined for convenience.
    Matrix transformation;
    transformation << t[0].x(), t[1].x(), t[2].x(),
                      t[0].y(), t[1].y(), t[2].y(),
                      t[0].z(), t[1].z(), t[2].z();

    FT area = CGAL::approximate_sqrt(t.squared_area());

    // skip zero measure primitives
    if(area == (FT)0.0)
      continue;

    // Find the 2nd order moment for the triangle wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = 2 * area * transformation * moment * transformation.transpose();

    // and add to covariance matrix
    covariance[0] += transformation(0,0);
    covariance[1] += transformation(1,0);
    covariance[2] += transformation(2,0);
    covariance[3] += transformation(1,1);
    covariance[4] += transformation(2,1);
    covariance[5] += transformation(2,2);

    mass += area;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += -mass * (c.x() * c.x());
  covariance[1] += -mass * (c.x() * c.y());
  covariance[2] += -mass * (c.z() * c.x());
  covariance[3] += -mass * (c.y() * c.y());
  covariance[4] += -mass * (c.z() * c.y());
  covariance[5] += -mass * (c.z() * c.z());

}

// assemble covariance matrix from a cuboid set
template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Eigen_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& ,                    // kernel
                             const typename K::Iso_cuboid_3*,// used for indirection
                             const CGAL::Dimension_tag<3>&,
                             const Eigen_diagonalize_traits<typename K::FT, 3>&)
{
  typedef typename K::FT          FT;
  typedef typename K::Iso_cuboid_3    Iso_cuboid;
  typedef typename Eigen::Matrix<FT, 3, 3> Matrix;

  // assemble covariance matrix as a semi-definite matrix.
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5
  // final combined covariance matrix for all cuboids and their combined mass
  FT mass = (FT)0.0;

  // assemble 2nd order moment about the origin.
  Matrix moment;
  moment << (FT)(1.0/3.0), (FT)(1.0/4.0), (FT)(1.0/4.0),
            (FT)(1.0/4.0), (FT)(1.0/3.0), (FT)(1.0/4.0),
            (FT)(1.0/4.0), (FT)(1.0/4.0), (FT)(1.0/3.0);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each cuboid, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Iso_cuboid& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT x0 = t.xmin();
    FT y0 = t.ymin();
    FT z0 = t.zmin();

    FT x1 = t.xmax();
    FT y1 = t.ymax();
    FT z1 = t.zmax();

    Matrix transformation;
    transformation << x1 - x0, 0      , 0      ,
                      0      , y1 - y0, 0      ,
                      0      , 0      , z1 - z0;

    FT volume = (x1-x0) * (y1-y0) * (z1-z0);

    // skip zero measure primitives
    if(volume == (FT)0.0)
      continue;
    CGAL_assertion(volume > 0.0);

    // Find the 2nd order moment for the cuboid wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = volume * transformation * moment * transformation.transpose();

    // Translate the 2nd order moment to the minimum corner (x0,y0,z0) of the cuboid.
    FT xav0 = (x1 - x0) / FT(2.0);
    FT yav0 = (y1 - y0) / FT(2.0);
    FT zav0 = (z1 - z0) / FT(2.0);

    // and add to covariance matrix
    covariance[0] += transformation(0,0) + volume * (2*x0*xav0 + x0*x0);
    covariance[1] += transformation(1,0) + volume * (xav0*y0 + yav0*x0 + x0*y0);
    covariance[2] += transformation(2,0) + volume * (x0*zav0 + xav0*z0 + x0*z0);
    covariance[3] += transformation(1,1) + volume * (2*y0*yav0 + y0*y0);
    covariance[4] += transformation(2,1) + volume * (yav0*z0 + y0*zav0 + z0*y0);
    covariance[5] += transformation(2,2) + volume * (2*zav0*z0 + z0*z0);

    mass += volume;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

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
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Eigen_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& ,                    // kernel
                             const typename K::Iso_cuboid_3*,// used for indirection
                             const CGAL::Dimension_tag<2>&,
                             const Eigen_diagonalize_traits<typename K::FT, 3>&)
{
#if 1
  typedef typename K::FT              FT;
  typedef typename K::Iso_cuboid_3    Iso_cuboid;
  typedef typename K::Triangle_3      Triangle;
  auto converter = [](const Iso_cuboid& c, int idx) -> Triangle
    {
      // Decomposition of 6 faces of the cuboid into 12 triangles
      static constexpr std::array<std::array<int, 3>, 12 > indices
      = {{ { 0, 1, 2 }, { 0, 2, 3 }, { 2, 3, 4 }, { 2, 4, 7 },
           { 3, 4, 5 }, { 3, 5, 0 }, { 4, 5, 6 }, { 4, 6, 7 },
           { 5, 6, 1 }, { 5, 1, 0 }, { 6, 7, 2 }, { 6, 2, 1 } }};
      return Triangle (c[indices[idx][0]], c[indices[idx][1]], c[indices[idx][2]]);
    };

  assemble_covariance_matrix_3
    (make_subiterator<Triangle, 12> (first, converter),
     make_subiterator<Triangle, 12> (beyond),
     covariance, c, K(), (Triangle*)nullptr, CGAL::Dimension_tag<2>(),
     Eigen_diagonalize_traits<FT, 3>());

#else
  // This variant uses the standard formulas but seems to be broken
  // (line/plane estimated appear to be wrong). In the absence of a
  // reliable fix so far, the above workaround applying PCA to a
  // decomposition of the cuboid into triangles is used.

  typedef typename K::FT FT;
  typedef typename K::Iso_cuboid_3 Iso_cuboid;
  typedef typename Eigen::Matrix<FT, 3, 3> Matrix;

  // assemble covariance matrix as a semi-definite matrix.
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5
  //Final combined covariance matrix for all cuboids and their combined mass
  FT mass = (FT)0.0;

  // assemble 2nd order moment about the origin.
  Matrix moment;
  moment << (FT)(7.0/3.0), (FT)1.5,       (FT)1.5,
            (FT)1.5,       (FT)(7.0/3.0), (FT)1.5,
            (FT)1.5,       (FT)1.5,       (FT)(7.0/3.0);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each cuboid, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Iso_cuboid& t = *it;

    // defined for convenience.
    FT x0 = t.xmin();
    FT y0 = t.ymin();
    FT z0 = t.zmin();

    FT x1 = t.xmax();
    FT y1 = t.ymax();
    FT z1 = t.zmax();

    Matrix transformation;
    transformation << x1 - x0, 0      , 0      ,
                      0      , y1 - y0, 0      ,
                      0      , 0      , z1 - z0;

    FT area = FT(2) * ((x1-x0)*(y1-y0) + (x1-x0)*(z1-z0) + (y1-y0)*(z1-z0));

    // skip zero measure primitives
    if(area == (FT)0.0)
      continue;
    CGAL_assertion(area > 0.0);

    // Find the 2nd order moment for the cuboid wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = area * transformation * moment * transformation.transpose();

    // Translate the 2nd order moment to the minimum corner (x0,y0,z0) of the cuboid.
    FT xav0 = (x1 - x0) / (2.0);
    FT yav0 = (y1 - y0) / (2.0);
    FT zav0 = (z1 - z0) / (2.0);

    // and add to covariance matrix
    covariance[0] += transformation(0,0) + area * (2*x0*xav0 + x0*x0);
    covariance[1] += transformation(1,0) + area * (xav0*y0 + yav0*x0 + x0*y0);
    covariance[2] += transformation(2,0) + area * (x0*zav0 + xav0*z0 + x0*z0);
    covariance[3] += transformation(1,1) + area * (2*y0*yav0 + y0*y0);
    covariance[4] += transformation(2,1) + area * (yav0*z0 + y0*zav0 + z0*y0);
    covariance[5] += transformation(2,2) + area * (2*zav0*z0 + z0*z0);

    mass += area;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += -mass * (c.x() * c.x());
  covariance[1] += -mass * (c.x() * c.y());
  covariance[2] += -mass * (c.z() * c.x());
  covariance[3] += -mass * (c.y() * c.y());
  covariance[4] += -mass * (c.z() * c.y());
  covariance[5] += -mass * (c.z() * c.z());
#endif
}

// assemble covariance matrix from a sphere set
template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Eigen_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K&,                     // kernel
                             const typename K::Sphere_3*,  // used for indirection
                             const CGAL::Dimension_tag<3>&,
                             const Eigen_diagonalize_traits<typename K::FT, 3>&)
{
  typedef typename K::FT          FT;
  typedef typename K::Sphere_3  Sphere;
  typedef typename Eigen::Matrix<FT, 3, 3> Matrix;

  // assemble covariance matrix as a semi-definite matrix.
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5
  //Final combined covariance matrix for all spheres and their combined mass
  FT mass = 0.0;

  // assemble 2nd order moment about the origin.
  Matrix moment;
  moment << FT(4.0/15.0), FT(0.0),      FT(0.0),
            FT(0.0),      FT(4.0/15.0), FT(0.0),
            FT(0.0),      FT(0.0),      FT(4.0/15.0);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each sphere, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Sphere& t = *it;

    // defined for convenience.
    FT radius = CGAL::approximate_sqrt(t.squared_radius());
    Matrix transformation;
    transformation << radius, 0.0, 0.0,
                      0.0, radius, 0.0,
                      0.0, 0.0, radius;
    FT volume =  radius * t.squared_radius();

    // skip zero measure primitives
    if(volume == (FT)0.0)
      continue;

    // Find the 2nd order moment for the sphere wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = volume * transformation * moment * transformation.transpose();

    volume *= FT(4.0 / 3.0);
    // Translate the 2nd order moment to the center of the sphere.
    FT x0 = t.center().x();
    FT y0 = t.center().y();
    FT z0 = t.center().z();

    // and add to covariance matrix
    covariance[0] += transformation(0,0) + volume * x0*x0;
    covariance[1] += transformation(1,0) + volume * x0*y0;
    covariance[2] += transformation(2,0) + volume * x0*z0;
    covariance[3] += transformation(1,1) + volume * y0*y0;
    covariance[4] += transformation(2,1) + volume * z0*y0;
    covariance[5] += transformation(2,2) + volume * z0*z0;

    mass += volume;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += -mass * (c.x() * c.x());
  covariance[1] += -mass * (c.x() * c.y());
  covariance[2] += -mass * (c.z() * c.x());
  covariance[3] += -mass * (c.y() * c.y());
  covariance[4] += -mass * (c.z() * c.y());
  covariance[5] += -mass * (c.z() * c.z());

}
// assemble covariance matrix from a sphere set
template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Eigen_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K&,                     // kernel
                             const typename K::Sphere_3*,  // used for indirection
                             const CGAL::Dimension_tag<2>&,
                             const Eigen_diagonalize_traits<typename K::FT, 3>&)
{
  typedef typename K::FT          FT;
  typedef typename K::Sphere_3  Sphere;
  typedef typename Eigen::Matrix<FT, 3, 3> Matrix;

  // assemble covariance matrix as a semi-definite matrix.
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5
  //Final combined covariance matrix for all spheres and their combined mass
  FT mass = 0.0;

  // assemble 2nd order moment about the origin.
  Matrix moment;
  moment << FT(4.0/3.0), FT(0.0),     FT(0.0),
            FT(0.0),     FT(4.0/3.0), FT(0.0),
            FT(0.0),     FT(0.0),     FT(4.0/3.0);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each sphere, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Sphere& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT radius = CGAL::approximate_sqrt(t.squared_radius());
    Matrix transformation;
    transformation << radius, 0.0,    0.0,
                      0.0,    radius, 0.0,
                      0.0,    0.0,    radius;
    FT area = t.squared_radius();

    // skip zero measure primitives
    if(area == (FT)0.0)
      continue;

    // Find the 2nd order moment for the sphere wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = area * transformation * moment * transformation.transpose();

    area *= FT(4.0);
    // Translate the 2nd order moment to the center of the sphere.
    FT x0 = t.center().x();
    FT y0 = t.center().y();
    FT z0 = t.center().z();

    // and add to covariance matrix
    covariance[0] += transformation(0,0) + area * x0*x0;
    covariance[1] += transformation(1,0) + area * x0*y0;
    covariance[2] += transformation(2,0) + area * x0*z0;
    covariance[3] += transformation(1,1) + area * y0*y0;
    covariance[4] += transformation(2,1) + area * z0*y0;
    covariance[5] += transformation(2,2) + area * z0*z0;

    mass += area;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += -mass * (c.x() * c.x());
  covariance[1] += -mass * (c.x() * c.y());
  covariance[2] += -mass * (c.z() * c.x());
  covariance[3] += -mass * (c.y() * c.y());
  covariance[4] += -mass * (c.z() * c.y());
  covariance[5] += -mass * (c.z() * c.z());

}

// assemble covariance matrix from a tetrahedron set
template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Eigen_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& ,                    // kernel
                             const typename K::Tetrahedron_3*,// used for indirection
                             const CGAL::Dimension_tag<3>&,
                             const Eigen_diagonalize_traits<typename K::FT, 3>&)
{
  typedef typename K::FT          FT;
  typedef typename K::Point_3     Point_3;
  typedef typename K::Vector_3    Vector_3;
  typedef typename K::Tetrahedron_3  Tetrahedron;
  typedef typename Eigen::Matrix<FT, 3, 3> Matrix;
  typedef typename Eigen::Matrix<FT, 3, 1> Vector;

  // assemble covariance matrix as a semi-definite matrix.
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5
  // assemble 2nd order moment about the origin.
  Matrix moment;
  moment << FT(1.0/60.0),  FT(1.0/120.0), FT(1.0/120.0),
            FT(1.0/120.0), FT(1.0/60.0),  FT(1.0/120.0),
            FT(1.0/120.0), FT(1.0/120.0), FT(1.0/60.0);

  Matrix accum; // zero by default
  accum << 0, 0, 0, 0, 0, 0, 0, 0, 0;
  for (InputIterator it = first;  it != beyond; it++)
  {
    const Tetrahedron& t = *it;

    // defined for convenience.
    FT x0 = t[0].x();
    FT y0 = t[0].y();
    FT z0 = t[0].z();

    Matrix transformation;
    transformation << t[1].x()-x0, t[2].x()-x0, t[3].x()-x0,
                      t[1].y()-y0, t[2].y()-y0, t[3].y()-y0,
                      t[1].z()-z0, t[2].z()-z0, t[3].z()-z0;
    FT volume = CGAL::abs(t.volume());

    // skip zero measure primitives
    if(volume == (FT)0.0)
      continue;

    // affine transform
    transformation = 6. * volume * transformation * moment * transformation.transpose();

    Vector_3 d = t[0] - c; // delta
    Vector vec_d;
    vec_d << d.x(), d.y(), d.z();

    Point_3 C = CGAL::centroid(t) - (t[0] - CGAL::ORIGIN); // careful, local centroid
    Vector vec_c;
    vec_c << C.x(), C.y(), C.z();

    Matrix M = vec_c * vec_d.transpose() + vec_d * vec_c.transpose() + vec_d * vec_d.transpose();

    accum += transformation + volume * M;
  }

  covariance[0] = accum(0,0);
  covariance[1] = accum(1,0);
  covariance[2] = accum(2,0);
  covariance[3] = accum(1,1);
  covariance[4] = accum(2,1);
  covariance[5] = accum(2,2);

}

// assemble covariance matrix from a segment set
template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Eigen_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& ,                    // kernel
                             const typename K::Segment_3*,// used for indirection
                             const CGAL::Dimension_tag<1>&,
                             const Eigen_diagonalize_traits<typename K::FT, 3>&)
{
  typedef typename K::FT          FT;
  typedef typename K::Segment_3  Segment;
  typedef typename Eigen::Matrix<FT, 3, 3> Matrix;

  // assemble covariance matrix as a semi-definite matrix.
  // Matrix numbering:
  // 0 1 2
  //   3 4
  //     5
  //Final combined covariance matrix for all segments and their combined mass
  FT mass = 0.0;

  // assemble 2nd order moment about the origin.
  Matrix moment;
  moment << FT(1.0/3.0), FT(0.5/3.0), FT(0.0),
            FT(0.5/3.0), FT(1.0/3.0), FT(0.0),
            FT(0.0),     FT(0.0),     FT(0.0);

  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    // Now for each segment, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Segment& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    Matrix transformation;
    transformation << t[0].x(), t[1].x(), 0.0,
                      t[0].y(), t[1].y(), 0.0,
                      t[0].z(), t[1].z(), 1.0;
    FT length = CGAL::approximate_sqrt(t.squared_length());

    // skip zero measure primitives
    if(length == (FT)0.0)
      continue;

    // Find the 2nd order moment for the segment wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = length * transformation * moment * transformation.transpose();

    // and add to covariance matrix
    covariance[0] += transformation(0,0);
    covariance[1] += transformation(1,0);
    covariance[2] += transformation(2,0);
    covariance[3] += transformation(1,1);
    covariance[4] += transformation(2,1);
    covariance[5] += transformation(2,2);

    mass += length;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  covariance[0] += -mass * (c.x() * c.x());
  covariance[1] += -mass * (c.x() * c.y());
  covariance[2] += -mass * (c.z() * c.x());
  covariance[3] += -mass * (c.y() * c.y());
  covariance[4] += -mass * (c.z() * c.y());
  covariance[5] += -mass * (c.z() * c.z());

}

// Variants using default
template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Default_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& k,                    // kernel
                             const typename K::Triangle_3* t,// used for indirection
                             const CGAL::Dimension_tag<2>& tag,
                             const Default_diagonalize_traits<typename K::FT, 3>&)
{
  assemble_covariance_matrix_3 (first, beyond, covariance, c, k, t, tag,
                                Eigen_diagonalize_traits<typename K::FT, 3>());
}

template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Default_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& k,                    // kernel
                             const typename K::Iso_cuboid_3* ic,// used for indirection
                             const CGAL::Dimension_tag<3>& tag,
                             const Default_diagonalize_traits<typename K::FT, 3>&)
{
  assemble_covariance_matrix_3 (first, beyond, covariance, c, k, ic, tag,
                                Eigen_diagonalize_traits<typename K::FT, 3>());
}

template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Default_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& k,                    // kernel
                             const typename K::Iso_cuboid_3* ic,// used for indirection
                             const CGAL::Dimension_tag<2>& tag,
                             const Default_diagonalize_traits<typename K::FT, 3>&)
{
  assemble_covariance_matrix_3 (first, beyond, covariance, c, k, ic, tag,
                                Eigen_diagonalize_traits<typename K::FT, 3>());
}

template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Default_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& k,                     // kernel
                             const typename K::Sphere_3* s,  // used for indirection
                             const CGAL::Dimension_tag<3>& tag,
                             const Default_diagonalize_traits<typename K::FT, 3>&)
{
  assemble_covariance_matrix_3 (first, beyond, covariance, c, k, s, tag,
                                Eigen_diagonalize_traits<typename K::FT, 3>());
}


template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Default_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& k,                     // kernel
                             const typename K::Sphere_3* s,  // used for indirection
                             const CGAL::Dimension_tag<2>& tag,
                             const Default_diagonalize_traits<typename K::FT, 3>&)
{
  assemble_covariance_matrix_3 (first, beyond, covariance, c, k, s, tag,
                                Eigen_diagonalize_traits<typename K::FT, 3>());
}

template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Default_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& k,                    // kernel
                             const typename K::Tetrahedron_3* t,// used for indirection
                             const CGAL::Dimension_tag<3>& tag,
                             const Default_diagonalize_traits<typename K::FT, 3>&)
{
  assemble_covariance_matrix_3 (first, beyond, covariance, c, k, t, tag,
                                Eigen_diagonalize_traits<typename K::FT, 3>());
}

template < typename InputIterator,
           typename K >
void
assemble_covariance_matrix_3(InputIterator first,
                             InputIterator beyond,
                             typename Default_diagonalize_traits<typename K::FT, 3>::Covariance_matrix& covariance, // covariance matrix
                             const typename K::Point_3& c, // centroid
                             const K& k,                    // kernel
                             const typename K::Segment_3* s,// used for indirection
                             const CGAL::Dimension_tag<1>& tag,
                             const Default_diagonalize_traits<typename K::FT, 3>&)
{
  assemble_covariance_matrix_3 (first, beyond, covariance, c, k, s, tag,
                                Eigen_diagonalize_traits<typename K::FT, 3>());
}




} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_UTIL_EIGEN_H
