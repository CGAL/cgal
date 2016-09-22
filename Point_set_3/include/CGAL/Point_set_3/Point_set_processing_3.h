// Copyright (c) 2016  GeometryFactory Sarl (France).
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
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_POINT_SET_3_POINT_SET_PROCESSING_3_H
#define CGAL_POINT_SET_3_POINT_SET_PROCESSING_3_H

#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/edge_aware_upsample_point_set.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/hierarchy_simplify_point_set.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/remove_outliers.h>
#include <CGAL/vcm_estimate_normals.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/read_ply_point_set_3.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/IO/write_off_points.h>
#include <CGAL/IO/write_ply_points.h>


namespace CGAL {

template <typename Concurrency_tag,
          typename PointSet>
double
bilateral_smooth_point_set(
  PointSet& point_set, ///< point set
  const unsigned int k,     ///< number of neighbors.
  double sharpness_angle)   ///< control sharpness(0-90)
{
  return CGAL::bilateral_smooth_point_set<Concurrency_tag>
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(),
     k, sharpness_angle);
}


template <typename Concurrency_tag,
	  typename PointSet>
double
compute_average_spacing(
  const PointSet& point_set, ///< point set
  unsigned int k) ///< number of neighbors.
{
  return CGAL::compute_average_spacing<Concurrency_tag>
    (point_set.begin(), point_set.end(), point_set.point_map(), k);
}

template <typename Concurrency_tag,
          typename PointSet>
void
edge_aware_upsample_point_set(
  PointSet& point_set, ///< point set
  double sharpness_angle = 30,  ///< control sharpness(0-90)
  double edge_sensitivity = 1,  ///< edge sensitivity(0-5)
  double neighbor_radius = -1, ///< initial size of neighbors.
  const std::size_t number_of_output_points = 1000)///< number of output points.     
{
  CGAL::edge_aware_upsample_point_set<Concurrency_tag>
    (point_set.begin(), point_set.end(),
     point_set.point_back_inserter(),
     point_set.point_map(), point_set.normal_map(),
     sharpness_angle, edge_sensitivity, neighbor_radius, number_of_output_points);
}

template <typename PointSet>
void grid_simplify_point_set(
  PointSet& point_set, ///< point set
  double epsilon) ///< tolerance value when merging 3D points.
{
  point_set.remove_from
    (CGAL::grid_simplify_point_set
     (point_set.begin(), point_set.end(), point_set.point_map(), epsilon));
}

template <typename PointSet>
void hierarchy_simplify_point_set(
  PointSet& point_set, ///< point set
  const unsigned int size = 10,
  const double var_max = 0.333)
{
  point_set.remove_from
    (CGAL::hierarchy_simplify_point_set
     (point_set.begin(), point_set.end(), point_set.point_map(), size, var_max,
      CGAL::Default_diagonalize_traits<double, 3>()));
}

  
template <typename Concurrency_tag,
	  typename PointSet>
void
jet_estimate_normals(
  PointSet& point_set, ///< point set
  unsigned int k, ///< number of neighbors.
  unsigned int degree_fitting = 2) ///< fitting degree
{
  point_set.add_normal_property();

  CGAL::jet_estimate_normals<Concurrency_tag>
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(),
     k, degree_fitting);
}

  
template <typename Concurrency_tag,
	  typename PointSet>
void
jet_smooth_point_set(
  PointSet& point_set, ///< point set
  unsigned int k, ///< number of neighbors.
  unsigned int degree_fitting = 2, ///< fitting degree
  unsigned int degree_monge = 2)  ///< Monge degree
{
  CGAL::jet_smooth_point_set<Concurrency_tag>
    (point_set.begin(), point_set.end(), point_set.point_map(),
     k, degree_fitting, degree_monge);

}

template <typename PointSet>
typename PointSet::iterator
mst_orient_normals(
  PointSet& point_set, ///< point set
  unsigned int k) ///< number of neighbors
{
  return CGAL::mst_orient_normals
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(), k);
}

template <typename Concurrency_tag,
	  typename PointSet>
void
pca_estimate_normals(
  PointSet& point_set, ///< point set
  unsigned int k) ///< number of neighbors.
{
  point_set.add_normal_property();

  CGAL::pca_estimate_normals<Concurrency_tag>
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(),
     k);
}

template <typename PointSet>
void random_simplify_point_set(
  PointSet& point_set, ///< point set
  double removed_percentage) ///< percentage of points to remove
{
  point_set.remove_from
    (CGAL::random_simplify_point_set
     (point_set.begin(), point_set.end(), point_set.point_map(), removed_percentage));
}


template <typename PointSet>
void remove_outliers(
  PointSet& point_set, ///< point set
  unsigned int k, ///< number of neighbors.
  double threshold_percent) ///< percentage of points to remove
{
  point_set.remove_from
    (CGAL::remove_outliers
     (point_set.begin(), point_set.end(), point_set.point_map(), k, threshold_percent));
}

template <typename PointSet>
void
vcm_estimate_normals(
  PointSet& point_set, ///< point set
  double offset_radius,
  double convolution_radius)
{
  point_set.add_normal_property();

  CGAL::vcm_estimate_normals
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(),
     offset_radius, convolution_radius);
}

template <typename PointSet>
void
vcm_estimate_normals(
  PointSet& point_set, ///< point set
  double offset_radius,
  unsigned int nb_neighbors_convolve)
{
  point_set.add_normal_property();

  CGAL::vcm_estimate_normals
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(),
     offset_radius, nb_neighbors_convolve);
}

template <typename Concurrency_tag,
          typename PointSet>
void
wlop_simplify_and_regularize_point_set(
  const PointSet& input_point_set, ///< input point set
  PointSet& output_point_set, ///< output point set
  const double select_percentage = 5,     ///< percentage of points to retain
  double neighbor_radius = -1,       ///< size of neighbors.
  const unsigned int max_iter_number = 35, ///< number of iterations.
  const bool require_uniform_sampling = false     ///< if needed to compute density 
                                      ///  to generate more rugularized result.                                 
) 
{
  CGAL::wlop_simplify_and_regularize_point_set
    (input_point_set.begin(), input_point_set.end(),
     output_point_set.point_back_inserter(),
     input_point_set.point_map(),
     neighbor_radius, max_iter_number, require_uniform_sampling);
}

template <typename PointSet>
bool
read_xyz_point_set(
  std::istream& stream, ///< input stream.
  PointSet& point_set) ///< point set
{
  point_set.add_normal_property();

  bool out = CGAL::read_xyz_points_and_normals
    (stream,
     point_set.index_back_inserter(),
     point_set.point_push_map(),
     point_set.normal_push_map());

  bool has_normals = false;
  for (typename PointSet::const_iterator it = point_set.begin();
       it != point_set.end(); ++ it)
    if (point_set.normal(it) != CGAL::NULL_VECTOR)
      {
        has_normals = true;
        break;
      }

  if (!has_normals)
    point_set.remove_normal_property();
  
  return out;
}

template <typename PointSet>
bool
read_off_point_set(
  std::istream& stream, ///< input stream.
  PointSet& point_set) ///< point set
{
  point_set.add_normal_property();

  bool out = CGAL::read_off_points_and_normals
    (stream,
     point_set.index_back_inserter(),
     point_set.point_push_map(),
     point_set.normal_push_map());

  bool has_normals = false;
  for (typename PointSet::const_iterator it = point_set.begin();
       it != point_set.end(); ++ it)
    if (point_set.normal(it) != CGAL::NULL_VECTOR)
      {
        has_normals = true;
        break;
      }

  if (!has_normals)
    point_set.remove_normal_property();

  return out;
}

template <typename PointSet>
bool
read_ply_point_set(
  std::istream& stream, ///< input stream.
  PointSet& point_set) ///< point set
{
  
  typedef typename PointSet::Point_type Point;
  typedef typename PointSet::Vector_type Vector;

  CGAL::Ply_interpreter_point_set_3<Point, Vector> interpreter (point_set);

  return CGAL::read_ply_custom_points
    (stream, interpreter,
     typename Kernel_traits<Point>::Kernel());
}

template <typename PointSet>
bool
write_xyz_point_set(
  std::ostream& stream, ///< output stream.
  const PointSet& point_set)  ///< point set
{
  if (point_set.has_normals())
    return CGAL::write_xyz_points_and_normals
      (stream, point_set.begin(), point_set.end(),
       point_set.point_map(), point_set.normal_map());
  
  return CGAL::write_xyz_points
  (stream, point_set.begin(), point_set.end(),
   point_set.point_map());
}

template <typename PointSet>
bool
write_off_point_set(
  std::ostream& stream, ///< output stream.
  const PointSet& point_set)  ///< point set
{
  if (point_set.has_normals())
    return CGAL::write_off_points_and_normals
      (stream, point_set.begin(), point_set.end(),
       point_set.point_map(), point_set.normal_map());
  
  return CGAL::write_off_points
  (stream, point_set.begin(), point_set.end(),
   point_set.point_map());
}

template <typename PointSet>
bool
write_ply_point_set(
  std::ostream& stream, ///< output stream.
  const PointSet& point_set)  ///< point set
{
  if (point_set.has_normals())
    return CGAL::write_ply_points_and_normals
      (stream, point_set.begin(), point_set.end(),
       point_set.point_map(), point_set.normal_map());
  
  return CGAL::write_ply_points
  (stream, point_set.begin(), point_set.end(),
   point_set.point_map());
}

} // namespace CGAL


#endif // CGAL_POINT_SET_3_POINT_SET_PROCESSING_3_H
