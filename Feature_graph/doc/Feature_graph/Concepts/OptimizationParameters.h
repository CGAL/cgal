/*!
\ingroup PkgFeatureGraphConcepts
\cgalConcept

The concept 'OptimizationParameters' describes the parameters
for the optimization step.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Optimization_parameters}
\cgalHasModelsEnd

*/

class OptimizationParameters {
public:


/// \name Types
/// @{

/*!
Natural number type.
*/
typedef unspecified_type Size;

/*!
Numerical type.
*/
typedef unspecified_type FT;

/// @}

/// \name Access Functions
/// @{

/*!
Returns the maximum number of iteration of the gradient descent.
*/
Size maximum_number_of_iteration() const;
/*!
Returns the step size at the first iteration of the gradient descent.
*/
FT start_step_size() const;
/*!
Returns the step size at the last iteration of the gradient descent.
*/
FT end_step_size() const;
/*!
Returns the minimum energy change to stop the gradient descent iterations.
*/
FT mininmum_energy_delta() const;
/*!
Returns the distance to collapse adjacent points in a line during the gradient descent.
*/
FT collapse_distance() const;
/*!
Returns the smoothing factor of the energy.
0 means no smoothing,
1 means that the energy will consider smoothing with the same weight
as the displacement toward the sharp features of the surface.
*/
FT smooth_factor() const;
/*!
Returns the distance to refine the normals of elements near the sharp features.
*/
FT refine_normal_distance() const;
/*!
Returns the distance to collect elements near the sharp features
to determine the adjacent planes.
*/
FT plane_detection_distance() const;
/*!
Returns a functor that estimates the normal on an element.

\tparam Normal_estimator a model of `NormalEstimator`.
*/
Normal_estimator normal_estimator() const;

/// @}

}; /* end OptimizationParameters */