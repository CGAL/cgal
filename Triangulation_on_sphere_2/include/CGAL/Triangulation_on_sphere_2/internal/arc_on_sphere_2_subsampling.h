// Copyright (c) 1997-2021 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sébastien Loriot
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_TOS2_INTERNAL_ARC_ON_SPHERE_SUBSAMPLING_H
#define CGAL_TOS2_INTERNAL_ARC_ON_SPHERE_SUBSAMPLING_H

#include <CGAL/license/Triangulation_on_sphere_2.h>

#include <CGAL/assertions.h>
#include <CGAL/Default.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <cmath>
#include <list>

namespace CGAL {
namespace Triangulations_on_sphere_2 {
namespace internal {

template <class Kernel,
          class Matrix_ = Default,
          class Col_ = Default,
          class EigenlessDefault = void>
double get_theta(typename Kernel::Point_3& pt,
                 typename Kernel::Vector_3& V1,
                 typename Kernel::Vector_3& V2,
                 typename Kernel::Vector_3& V3)
{
  typedef typename Kernel::FT                                        FT;

  typedef typename Default::Get<Matrix_,
#ifdef CGAL_EIGEN3_ENABLED
                                Eigen::Matrix<FT, 3, 3, Eigen::DontAlign>
#else
                                EigenlessDefault
#endif
                                >::type                              Matrix;

  typedef typename Default::Get<Col_,
#ifdef CGAL_EIGEN3_ENABLED
                                Eigen::Matrix<FT, 3, 1>
#else
                                EigenlessDefault
#endif
                                >::type                              Col;

  static_assert(!(std::is_same<Matrix, EigenlessDefault>::value),
                            "Eigen is required to perform arc subsampling!");

  auto V1c = V1.cartesian_begin(), V2c = V2.cartesian_begin(), V3c = V3.cartesian_begin();

  Matrix M;
  M(0, 0) = *V1c++; M(0, 1) = *V2c++; M(0, 2) = *V3c++;
  M(1, 0) = *V1c++; M(1, 1) = *V2c++; M(1, 2) = *V3c++;
  M(2, 0) = *V1c++; M(2, 1) = *V2c++; M(2, 2) = *V3c++;

  Col pt_in_vect;
  auto ptc = pt.cartesian_begin();
  pt_in_vect(0) = *ptc++; pt_in_vect(1) = *ptc++; pt_in_vect(2) = *ptc++;

  Matrix inverse = M.inverse();
  FT X = inverse.row(0) * pt_in_vect;
  FT Y = inverse.row(1) * pt_in_vect;

  double angle = atan2(Y*sqrt(V2.squared_length()), X*sqrt(V1.squared_length()));

  return angle;
}

template <class Kernel>
inline typename Kernel::Point_3 compute_point(const typename Kernel::Point_3& center,
                                              double radius,
                                              double theta,
                                              const typename Kernel::Vector_3& b1,
                                              const typename Kernel::Vector_3& b2)
{
  return center +
         radius * cos(theta) * b1 / sqrt(b1.squared_length()) +
         radius * sin(theta) * b2 / sqrt(b2.squared_length());
}

template <class Kernel, class Output_iterator>
void subsample_arc_on_sphere_2(const typename Kernel::Circle_3& circle,
                              double source,
                              double target,
                              const typename Kernel::Vector_3& b1,
                              const typename Kernel::Vector_3& b2,
                              Output_iterator out_pts,
                              double min_edge_size)
{
  if(source > target)
    target += 2*CGAL_PI;

  CGAL_assertion(target > source);

  const double radius = sqrt(circle.squared_radius());
  const double edge_len = (target - source) * radius;
  const int nb_of_segments = static_cast<int>(floor(edge_len / min_edge_size));

  *out_pts++ = compute_point<Kernel>(circle.center(), radius, source, b1, b2);
  const double step_size = (target - source) / static_cast<double>(nb_of_segments);
  double current_theta = source;
  for(int i=0; i<nb_of_segments-1; ++i)
  {
    current_theta += step_size;
    CGAL_assertion(current_theta <= target);
    *out_pts++ = compute_point<Kernel>(circle.center(), radius, current_theta, b1, b2);
  }
  *out_pts++ = compute_point<Kernel>(circle.center(), radius, target, b1, b2);
}

//Subsample from source to target seen ccw from the side of plane pointed by its orthogonal_vector()
template <class Kernel,class Output_iterator>
void subsample_circle_3(const typename Kernel::Circle_3& circle,
                        Output_iterator out_pts,
                        double min_edge_size)
{
  typename Kernel::Vector_3 b1 = circle.supporting_plane().base1();
  typename Kernel::Vector_3 b2 = circle.supporting_plane().base2();

  subsample_arc_on_sphere_2<Kernel>(circle, 0, 2*CGAL_PI, b1, b2, out_pts, min_edge_size);
}

//Subsample from source to target seen ccw from the side of plane pointed by its orthogonal_vector()
template <class Kernel,class Output_iterator>
void subsample_arc_on_sphere_2(const typename Kernel::Circle_3& circle,
                               const typename Kernel::Plane_3& plane,
                               const typename Kernel::Circular_arc_point_3& source,
                               const typename Kernel::Circular_arc_point_3& target,
                               Output_iterator out_pts,
                               double min_edge_size)
{
  typename Kernel::Vector_3 b1 = plane.base1();
  typename Kernel::Vector_3 b2 = plane.base2();
  typename Kernel::Vector_3 b3 = plane.orthogonal_vector();

  typename Kernel::Point_3 tmp = typename Kernel::Point_3(source.x() - circle.center().x(),
                                                          source.y() - circle.center().y(),
                                                          source.z() - circle.center().z());
  const double theta_source = get_theta<Kernel>(tmp, b1, b2, b3);

  tmp = typename Kernel::Point_3(target.x() - circle.center().x(),
                                 target.y() - circle.center().y(),
                                 target.z() - circle.center().z());
  const double theta_target = get_theta<Kernel>(tmp, b1, b2, b3);

  subsample_arc_on_sphere_2<Kernel>(circle, theta_source, theta_target, b1, b2, out_pts, min_edge_size);
}

template <class Kernel, class ArcOnSphere, class Output_iterator>
void subsample_arc_on_sphere_2(const ArcOnSphere& arc,
                               Output_iterator out_pts,
                               double min_edge_size = 0.01)
{
  subsample_arc_on_sphere_2<Kernel>(arc.supporting_circle(),
                                    arc.supporting_circle().supporting_plane(),
                                    arc.source(),
                                    arc.target(),
                                    out_pts,
                                    min_edge_size);
}

} // namespace internal
} // namespace Triangulations_on_sphere_2
} // namespace CGAL

#endif // CGAL_TOS2_INTERNAL_ARC_ON_SPHERE_SUBSAMPLING_H
