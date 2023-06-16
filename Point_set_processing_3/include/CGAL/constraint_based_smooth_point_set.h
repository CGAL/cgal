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
#include <CGAL/Point_set_processing_3/internal/Callback_wrapper.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

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

  // double neighbour_threshold = choose_parameter(get_parameter(np, internal_np::neighbour_threshold), 1.);
  // double normal_threshold = choose_parameter(get_parameter(np, internal_np::normal_threshold), 0.5);
  // double damping_factor = choose_parameter(get_parameter(np, internal_np::damping_factor), 3.);
  // double eigenvalue_threshold = choose_parameter(get_parameter(np, internal_np::eigenvalue_threshold), 1.);
  // double update_threshold = choose_parameter(get_parameter(np, internal_np::update_threshold), 2.);
  const std::function<bool(double)>& callback = choose_parameter(get_parameter(np, internal_np::callback),
                                                                 std::function<bool(double)>());

  double neighbour_threshold = 1.;
  double normal_threshold = 0.5;
  double damping_factor = 3.;
  double eigenvalue_threshold = 1.;
  double update_threshold = 2.;

  CGAL_precondition(points.begin() != points.end());

  // types for neighbors search structure
  typedef Point_set_processing_3::internal::Neighbor_query<Kernel, PointRange&, PointMap> Neighbor_query;

  PointMap point_map = NP_helper::get_point_map(points, np);
  NormalMap normal_map = NP_helper::get_normal_map(points, np);

  std::size_t nb_points = points.size();



  std::cout << "fn complete" << std::endl;

}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_CONSTRAINT_BASED_SMOOTH_POINT_SET_H
