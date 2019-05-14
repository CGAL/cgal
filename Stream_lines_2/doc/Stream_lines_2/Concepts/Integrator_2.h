
/*!
\ingroup PkgPlacementOfStreamlines2Concepts
\cgalConcept

The concept `Integrator_2` describes the set of requirements for 
the second template parameter of the class 
`CGAL::Stream_lines_2<VectorField_2,Integrator_2>`. This concept 
provides the operation that integrates a new point from a given point 
with a predefined step, and according to a specified vector. 

\cgalHasModel `CGAL::Euler_integrator_2<VectorField_2>`
\cgalHasModel `CGAL::Runge_kutta_integrator_2<VectorField_2>` 

*/

class Integrator_2 {
public:

/// \name Types 
/// @{

/*!
The scalar type. 
*/ 
typedef unspecified_type FT; 

/*!
The point type. 
*/ 
typedef unspecified_type Point_2; 

/*!
The vector type. 
*/ 
typedef unspecified_type Vector_2; 

/*!
The vector field type. 
*/ 
typedef unspecified_type Vector_field_2; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor.
*/ 
Integrator_2(); 

/// @} 

/// \name Operations 
/// The following operations return the newly integrated point.
/// @{

/*!
Returns the new position from the actual position defined by `p`, according to the vector given by `vector_field_2` at `p`. 
\pre `vector_field_2.is_in_domain(p)` must be `true`. 
*/ 
Point_2 operator()(Point_2 p, Vector_field_2 vector_field_2); 

/*!
\brief As above.\ The integration step is defined by `integration_step`. 
\pre `vector_field_2.is_in_domain(p)` must be `true`. 
*/ 
Point_2 operator()(Point_2 p, Vector_field_2 vector_field_2, FT integration_step); 

/*!
\brief As above.\ In addition, this function integrates forward if `direction` is true, and backward if it is false. 
\pre `vector_field_2.is_in_domain(p)` must be `true`. 
*/ 
Point_2 operator()(Point_2 p, Vector_field_2 vector_field_2, FT integration_step, bool direction); 

/// @}

}; /* end Integrator_2 */

