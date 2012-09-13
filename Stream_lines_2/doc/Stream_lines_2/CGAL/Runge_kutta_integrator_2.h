namespace CGAL {

/*!
\ingroup PkgPlacementOfStreamlines2

The template parameter `VectorField_2` has to 
be instantiated by a model of the concept `VectorField_2`. This class implements the second order Runge-Kutta integrator. 

\models ::Integrator_2 

\sa `Euler_integrator_2<VectorField_2>` 

*/
template< typename VectorField_2 >
class Runge_kutta_integrator_2 {
public:

/// \name Creation 
/// @{

/*! 
Creates a Runge-Kutta second order integrator class `rkinteg` with `integration_step` as integration step. 
*/ 
Runge_kutta_integrator_2(const FT & integration_step); 

/// @}

}; /* end Runge_kutta_integrator_2 */
} /* end namespace CGAL */
