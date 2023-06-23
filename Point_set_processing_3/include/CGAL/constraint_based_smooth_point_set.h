// Copyright (c) 2013-06  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : 

#ifndef CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H
#define CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Point_set_processing_3/internal/Neighbor_query.h>
// #include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>
#include <CGAL/for_each.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/iterator/zip_iterator.hpp>

#include <CGAL/Eigen_vector.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Kernel/global_functions.h>

namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL

namespace internal {

// very high level for now


// construct_nvt(point, neighbours)
// constructs nvt for a point and its neighbours
// return nvt

// decompose_nvt(nvt)
// get eigenvalues/vectors and perform BEO
// use diagonalize_selfadjoint_covariance_matrix()
// return "optimized" nvt

// nvt_normal_denoising(points, optimized nvts)
// calculates new better normals
// modifies points?

// calculate_covariance_matrix(points, neighbors)
// returns covariance matrix

// feature_detection(point, covar. matrix)
// classifies points based on covar (or based on something)
// returns enum?, maybe combine with next fn

// update_position(point, classification)
// implements the formulas to update position
// modifies point

template <typename Kernel, typename PointRange,
          typename PointMap, typename VectorMap>
Eigen::MatrixXd construct_nvt(
  const typename PointRange::iterator::value_type& vt,
  PointMap point_map,
  VectorMap normal_map,
  const std::vector<typename PointRange::iterator>& neighbor_pwns,
  typename Kernel::FT normal_threshold)
{
  // basic geometric types
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Point_3 Point;

  unsigned int w = 0; // cumulative weight
  Eigen::MatrixXd nvt(3, 3);

  const Point& p = get(point_map, vt);
  const Vector& n = get(normal_map, vt);

  for (typename PointRange::iterator it : neighbor_pwns)
  {
    const Point& np = get(point_map, *it);
    const Vector& nn = get(normal_map, *it);

    FT angle_difference = CGAL::approximate_angle(n, nn);
    if(angle_difference >= normal_threshold) {
      w += 1;
      Eigen::Vector3d vn(n.x(), n.y(), n.z());
      Eigen::Vector3d vnn(nn.x(), nn.y(), nn.z());

      nvt += vn * vnn.transpose();
    }
  }

  std::cout << nvt << std::endl;

  return nvt;
}


} /* namespace internal */

/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

template <typename ConcurrencyTag,
          typename PointRange,
          typename NamedParameters = parameters::Default_named_parameters
>
void
constraint_based_smooth_point_set(
  PointRange& points,
  const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // basic geometric types
  typedef typename PointRange::iterator iterator;
  typedef typename iterator::value_type value_type;
  typedef Point_set_processing_3_np_helper<PointRange, NamedParameters> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;
  typedef typename NP_helper::Geom_traits Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;

  CGAL_assertion_msg(NP_helper::has_normal_map(points, np), "Error: no normal map");

  typedef typename Kernel::FT FT;

  FT neighbor_radius = 1.;
  FT normal_threshold = 90.;
  FT damping_factor = 3.;
  FT eigenvalue_threshold = 1.;
  FT update_threshold = 2.;

  CGAL_precondition(points.begin() != points.end());

  // types for neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  PointMap point_map = NP_helper::get_point_map(points, np);
  NormalMap normal_map = NP_helper::get_normal_map(points, np);

  std::size_t nb_points = points.size();

  // initiate a KD-tree search for points
  Neighbor_query neighbor_query (points, point_map);

  // compute all neighbors
  typedef std::vector<iterator> iterators;
  std::vector<iterators> pwns_neighbors;
  pwns_neighbors.resize(nb_points);

  typedef boost::zip_iterator<boost::tuple<iterator, typename std::vector<iterators>::iterator> > Zip_iterator;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple (points.begin(), pwns_neighbors.begin())),
                      boost::make_zip_iterator (boost::make_tuple (points.end(), pwns_neighbors.end()))),
    [&](const typename Zip_iterator::reference& t)
    {
      neighbor_query.get_iterators (get(point_map, get<0>(t)), 0, neighbor_radius,
                                    std::back_inserter (get<1>(t)));
      return true;
    });

  std::vector<Eigen::MatrixXd> pwns_nvts(nb_points);

  typedef boost::zip_iterator
    <boost::tuple<iterator,
                  typename std::vector<iterators>::iterator,
                  typename std::vector<Eigen::MatrixXd>::iterator> > Zip_iterator_2;

  CGAL::for_each<ConcurrencyTag>
    (CGAL::make_range (boost::make_zip_iterator (boost::make_tuple
                                                (points.begin(), pwns_neighbors.begin(), pwns_nvts.begin())),
                      boost::make_zip_iterator (boost::make_tuple
                                                (points.end(), pwns_neighbors.end(), pwns_nvts.end()))),
    [&](const typename Zip_iterator_2::reference& t)
    {
      get<2>(t) = internal::construct_nvt<Kernel, PointRange>
          (get<0>(t),
           point_map, normal_map,
           get<1>(t),
           normal_threshold);
      return true;
    });
  

    



  std::cout << "fn complete" << std::endl;

}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H
