// Copyright (c) 2026 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ange Clement

#ifndef CGAL_OPTIMIZATION_PARAMETERS_H
#define CGAL_OPTIMIZATION_PARAMETERS_H

namespace CGAL
{

/*!
*
* \ingroup PkgFeatureGraphParameter
*
* The class `Optimization_parameters` describes the parameters for the optimization step.
* It contains constant values and a normal estimator.
*
* \tparam Normal_estimator a model of `NormalEstimator`.
*
* \sa `CGAL::Optimization_parameters_on_image`
* \sa `CGAL::Optimization_parameters_on_surface`
* \sa `CGAL::Detect_sharp_features_on_labeled_image`
* \sa `CGAL::Detect_sharp_features_on_surface`
*
*/

template <typename NormalEstimator>
class Optimization_parameters
{
public:

  /// \name Types
  /// @{

  /*!
  * Natural number type.
  */
  typedef std::size_t Size;

  /*!
  * Numerical type.
  */
  typedef double FT;

  /*!
  * Type of the functor that evaluates normals.
  */
  typedef NormalEstimator Normal_estimator;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  * constructs the parameters used to optimize the feature graph placement on a suface.
  *
  * \param maximum_iteration the maximum number of iteration of the gradient descent.
  * \param start_step_size the step size at the first iteration of the gradient descent.
  * \param end_step_size the step size at the last iteration of the gradient descent.
  * \param min_energy_delta the minimum energy change to stop the gradient descent iterations.
  * \param collapse_distance the distance to collapse adjacent points in a line during the gradient descent.
  * \param smooth_factor the smoothing factor of the energy.
  * 0 means no smoothing,
  * 1 means that the energy will consider smoothing with the same weight
  * as the displacement toward the sharp features of the surface.
  * \param refine_normal_distance the distance to refine the normals of elements near the sharp features.
  * \param plane_detection_distance the distance to collect elements near the sharp features
  * to determine the adjacent planes.
  * \param normal_estimator a functor to evaluate the normals on elements.
  */
  Optimization_parameters(
    const Size& maximum_iteration,
    const FT& start_step_size,
    const FT& end_step_size,
    const FT& min_energy_delta,
    const FT& collapse_distance,
    const FT& smooth_factor,
    const FT& refine_normal_distance,
    const FT& plane_detection_distance,
    const Normal_estimator& normal_estimator);

  /// @}

  /// \name Access Functions
  /// @{

  /*!
  * returns the maximum number of iteration of the gradient descent.
  */
  Size maximum_number_of_iteration() const;
  /*!
  * returns the step size at the first iteration of the gradient descent.
  */
  FT start_step_size() const;
  /*!
  * returns the step size at the last iteration of the gradient descent.
  */
  FT end_step_size() const;
  /*!
  * returns the minimum energy change to stop the gradient descent iterations.
  */
  FT mininmum_energy_delta() const;
  /*!
  * returns the distance to collapse adjacent points in a line during the gradient descent.
  */
  FT collapse_distance() const;
  /*!
  * returns the smoothing factor of the energy.
  * 0 means no smoothing,
  * 1 means that the energy will consider smoothing with the same weight
  * as the displacement toward the sharp features of the surface.
  */
  FT smooth_factor() const;
  /*!
  * returns the distance to refine the normals of elements near the sharp features.
  */
  FT refine_normal_distance() const;
  /*!
  * returns the distance to collect elements near the sharp features
  * to determine the adjacent planes.
  */
  FT plane_detection_distance() const;
  /*!
  * returns a functor that estimates the normal on an element.
  */
  Normal_estimator normal_estimator() const;

  /// @}

};

/*!
*
* \ingroup PkgFeatureGraphParameter
*
* The class `Optimization_parameters_on_image` describes the parameters for the optimization step
* with default values adapted for image inputs.
*
* \tparam Normal_estimator a model of `NormalEstimator`.
*
* \sa `CGAL::Optimization_parameters_on_surface`
* \sa `CGAL::Detect_sharp_features_on_labeled_image`
* \sa `CGAL::Detect_sharp_features_on_surface`
*
*/

template <typename NormalEstimator = CGAL::Normal_estimator::AmbrosioTortorelli_on_image>
class Optimization_parameters_on_image :
public Optimization_parameters<NormalEstimator>
{
private:
  typedef Optimization_parameters<NormalEstimator> Base;

public:

  /// \name Types
  /// @{

  /*!
  * Natural number type.
  */
  typedef typename Base::Size Size;

  /*!
  * Numerical type.
  */
  typedef typename Base::FT FT;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  * constructs the parameters used to optimize the feature graph placement on the suface of an image.
  *
  * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.
  *           The distances are expressed in terms of the longest voxel edge length.
  *
  * \cgalNamedParamsBegin
  *   \cgalParamSectionBegin{maximum_iteration}
  *     \cgalParamDescription{the maximum number of iteration of the gradient descent.}
  *     \cgalParamDefault{`Size(20)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{start_step_size}
  *     \cgalParamDescription{the step size at the first iteration of the gradient descent.}
  *     \cgalParamDefault{`FT(1.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{end_step_size}
  *     \cgalParamDescription{the step size at the last iteration of the gradient descent.}
  *     \cgalParamDefault{`FT(0.125)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{min_energy_delta}
  *     \cgalParamDescription{the minimum energy change to stop the gradient descent iterations.}
  *     \cgalParamDefault{`FT(1.e-3)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{collapse_distance}
  *     \cgalParamDescription{the distance to collapse adjacent points in a line during the gradient descent.}
  *     \cgalParamDefault{`FT(0.5)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{smooth_factor}
  *     \cgalParamDescription{the smoothing factor of the energy.
  *                           0 means no smoothing,
  *                           1 means that the energy will consider smoothing with the same weight
  *                           as the displacement toward the sharp features of the surface.}
  *     \cgalParamDefault{`FT(1.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{refine_normal_distance}
  *     \cgalParamDescription{the distance to refine the normals of elements near the sharp features.}
  *     \cgalParamDefault{`FT(4.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{plane_detection_distance}
  *     \cgalParamDescription{the distance to collect elements near the sharp features
  *                           to determine the adjacent planes.}
  *     \cgalParamDefault{`FT(4.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{normal_estimator}
  *     \cgalParamDescription{a functor to evaluate the normals on elements.}
  *     \cgalParamDefault{`CGAL::Feature_graph::AmbrosioTortorelli_on_image(image).get_normal_functor()`}
  *   \cgalParamSectionEnd
  * \cgalNamedParamsEnd
  *
  */
  template <typename CGAL_NP_TEMPLATE_PARAMETERS>
  Optimization_parameters_on_image(
    const CGAL_NP_CLASS& np = parameters::default_values());

  /// @}
};

/*!
*
* \ingroup PkgFeatureGraphParameter
*
* The class `Optimization_parameters_on_image` describes the parameters for the optimization step
* with default values adapted for surface inputs.
*
* \tparam Normal_estimator a model of `NormalEstimator`.
*
* \sa `CGAL::Optimization_parameters_on_image`
* \sa `CGAL::Detect_sharp_features_on_labeled_image`
* \sa `CGAL::Detect_sharp_features_on_surface`
*
*/

template <typename NormalEstimator = CGAL::Normal_estimator::Normal_estimator_on_surface>
class Optimization_parameters_on_surface :
public Optimization_parameters<NormalEstimator>
{
private:
  typedef Optimization_parameters<NormalEstimator> Base;

public:

  /// \name Types
  /// @{

  /*!
  * Natural number type.
  */
  typedef typename Base::Size Size;

  /*!
  * Numerical type.
  */
  typedef typename Base::FT FT;

  /// @}

  /// \name Constructor
  /// @{

  /*!
  * constructs the parameters used to optimize the feature graph placement on a suface.
  *
  * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.
  *
  * \cgalNamedParamsBegin
  *   \cgalParamSectionBegin{maximum_iteration}
  *     \cgalParamDescription{the maximum number of iteration of the gradient descent.}
  *     \cgalParamDefault{`Size(20)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{start_step_size}
  *     \cgalParamDescription{the step size at the first iteration of the gradient descent.}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{end_step_size}
  *     \cgalParamDescription{the step size at the last iteration of the gradient descent.}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{min_energy_delta}
  *     \cgalParamDescription{the minimum energy change to stop the gradient descent iterations.}
  *     \cgalParamDefault{`FT(1.e-3)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{collapse_distance}
  *     \cgalParamDescription{the distance to collapse adjacent points in a line during the gradient descent.}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{smooth_factor}
  *     \cgalParamDescription{the smoothing factor of the energy.
  *                           0 means no smoothing,
  *                           1 means that the energy will consider smoothing with the same weight
  *                           as the displacement toward the sharp features of the surface.}
  *     \cgalParamDefault{`FT(1.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{refine_normal_distance}
  *     \cgalParamDescription{the distance to refine the normals of elements near the sharp features.}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{plane_detection_distance}
  *     \cgalParamDescription{the distance to collect elements near the sharp features
  *                           to determine the adjacent planes.}
  *     \cgalParamDefault{`FT(0.0)`}
  *   \cgalParamSectionEnd
  *   \cgalParamSectionBegin{normal_estimator}
  *     \cgalParamDescription{a functor to evaluate the normals on elements.}
  *     \cgalParamDefault{`CGAL::Feature_graph::Normal_estimator::Normal_estimator_on_surface(surface)`}
  *   \cgalParamSectionEnd
  * \cgalNamedParamsEnd
  *
  */
  template <typename CGAL_NP_TEMPLATE_PARAMETERS>
  Optimization_parameters_on_surface(
    const CGAL_NP_CLASS& np = parameters::default_values());

  /// @}
};

} /* namespace CGAL */

#endif // CGAL_OPTIMIZATION_PARAMETERS_H