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

#include <CGAL/license/Point_set_3.h>


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


namespace CGAL {

/*!
  \ingroup PkgPointSet3PointSetProcessing3
 */
template <typename Concurrency_tag,
          typename Point,
          typename Vector>
double
bilateral_smooth_point_set(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  const unsigned int k,     ///< number of neighbors.
  double sharpness_angle)   ///< control sharpness(0-90)
{
  return CGAL::bilateral_smooth_point_set<Concurrency_tag>
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(),
     k, sharpness_angle);
}


/*!
  \ingroup PkgPointSet3PointSetProcessing3
 */
template <typename Concurrency_tag,
	  typename Point, typename Vector>
double
compute_average_spacing(
  const CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  unsigned int k) ///< number of neighbors.
{
  return CGAL::compute_average_spacing<Concurrency_tag>
    (point_set.begin(), point_set.end(), point_set.point_map(), k);
}

/*!
  \ingroup PkgPointSet3PointSetProcessing3
 */
template <typename Concurrency_tag,
          typename Point, typename Vector>
void
edge_aware_upsample_point_set(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
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

/*!
  \ingroup PkgPointSet3PointSetProcessing3

  \note No iterator is returned, points simplified are directly
  removed from the point set.
 */
template <typename Point, typename Vector>
void grid_simplify_point_set(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  double epsilon) ///< tolerance value when merging 3D points.
{
  point_set.remove_from
    (CGAL::grid_simplify_point_set
     (point_set.begin(), point_set.end(), point_set.point_map(), epsilon));
}

/*!
  \ingroup PkgPointSet3PointSetProcessing3

  \note No iterator is returned, points simplified are directly
  removed from the point set.
 */
template <typename Point, typename Vector>
void hierarchy_simplify_point_set(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  const unsigned int size = 10, ///< maximum cluster size
  const double var_max = 0.333) ///< maximal surface variation
{
  point_set.remove_from
    (CGAL::hierarchy_simplify_point_set
     (point_set.begin(), point_set.end(), point_set.point_map(), size, var_max,
      CGAL::Default_diagonalize_traits<double, 3>()));
}

  
/*!
  \ingroup PkgPointSet3PointSetProcessing3

  \note This function adds a normal map to the point set.
*/
template <typename Concurrency_tag,
	  typename Point, typename Vector>
void
jet_estimate_normals(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  unsigned int k, ///< number of neighbors.
  unsigned int degree_fitting = 2) ///< fitting degree
{
  point_set.add_normal_map();

  CGAL::jet_estimate_normals<Concurrency_tag>
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(),
     k, degree_fitting);
}

  
/*!
  \ingroup PkgPointSet3PointSetProcessing3
*/
template <typename Concurrency_tag,
	  typename Point, typename Vector>
void
jet_smooth_point_set(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  unsigned int k, ///< number of neighbors.
  unsigned int degree_fitting = 2, ///< fitting degree
  unsigned int degree_monge = 2)  ///< Monge degree
{
  CGAL::jet_smooth_point_set<Concurrency_tag>
    (point_set.begin(), point_set.end(), point_set.point_map(),
     k, degree_fitting, degree_monge);

}

/*!
  \ingroup PkgPointSet3PointSetProcessing3
*/
template <typename Point, typename Vector>
typename CGAL::Point_set_3<Point, Vector>::iterator
mst_orient_normals(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  unsigned int k) ///< number of neighbors
{
  return CGAL::mst_orient_normals
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(), k);
}

/*!
  \ingroup PkgPointSet3PointSetProcessing3

  \note This function adds a normal map to the point set.
*/
template <typename Concurrency_tag,
	  typename Point, typename Vector>
void
pca_estimate_normals(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  unsigned int k) ///< number of neighbors.
{
  point_set.add_normal_map();

  CGAL::pca_estimate_normals<Concurrency_tag>
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(),
     k);
}

/*!
  \ingroup PkgPointSet3PointSetProcessing3

  \note No iterator is returned, points simplified are directly
  removed from the point set.
 */
template <typename Point, typename Vector>
void random_simplify_point_set(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  double removed_percentage) ///< percentage of points to remove
{
  point_set.remove_from
    (CGAL::random_simplify_point_set
     (point_set.begin(), point_set.end(), point_set.point_map(), removed_percentage));
}


/*!
  \ingroup PkgPointSet3PointSetProcessing3

  \note No iterator is returned, points simplified are directly
  removed from the point set.
 */
template <typename Point, typename Vector>
void remove_outliers(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  unsigned int k, ///< number of neighbors.
  double threshold_percent) ///< percentage of points to remove
{
  point_set.remove_from
    (CGAL::remove_outliers
     (point_set.begin(), point_set.end(), point_set.point_map(), k, threshold_percent));
}

/*!
  \ingroup PkgPointSet3PointSetProcessing3

  \note This function adds a normal map to the point set.
*/
template <typename Point, typename Vector>
void
vcm_estimate_normals(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  double offset_radius, ///< offset radius.
  double convolution_radius) ///< convolution radius.
{
  point_set.add_normal_map();

  CGAL::vcm_estimate_normals
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(),
     offset_radius, convolution_radius);
}

/*!
  \ingroup PkgPointSet3PointSetProcessing3

  \note This function adds a normal map to the point set.
*/
template <typename Point, typename Vector>
void
vcm_estimate_normals(
  CGAL::Point_set_3<Point, Vector>& point_set, ///< point set
  double offset_radius, ///< offset radius.
  unsigned int nb_neighbors_convolve) ///< number of neighbors used during the convolution.
{
  point_set.add_normal_map();

  CGAL::vcm_estimate_normals
    (point_set.begin(), point_set.end(),
     point_set.point_map(), point_set.normal_map(),
     offset_radius, nb_neighbors_convolve);
}

/*!
  \ingroup PkgPointSet3PointSetProcessing3
*/
template <typename Concurrency_tag,
          typename Point, typename Vector>
void
wlop_simplify_and_regularize_point_set(
  const CGAL::Point_set_3<Point, Vector>& input_point_set, ///< input point set
  CGAL::Point_set_3<Point, Vector>& output_point_set, ///< output point set
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
     select_percentage,
     neighbor_radius,
     max_iter_number,
     require_uniform_sampling);
}


} // namespace CGAL


#endif // CGAL_POINT_SET_3_POINT_SET_PROCESSING_3_H
