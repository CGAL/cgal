namespace CGAL {

/*!
\ingroup PkgPlacementOfStreamlines2

This class implements the first order Euler integrator.

\tparam VectorField_2 has to be instantiated by a model of the concept `VectorField_2`. 

\cgalModels `Integrator_2`

\sa `Runge_kutta_integrator_2<VectorField_2>` 

*/
template< typename VectorField_2 >
class Euler_integrator_2 {
public:

/// \name Creation 
/// @{

/*!
Creates an Euler integrator with `integration_step` as integration step. 
*/ 
Euler_integrator_2(const FT & integration_step); 

/// @}

}; /* end Euler_integrator_2 */
} /* end namespace CGAL */
