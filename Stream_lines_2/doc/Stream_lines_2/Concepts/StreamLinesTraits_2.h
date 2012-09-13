
/*!
\ingroup PkgPlacementOfStreamlines2Concepts
\cgalconcept

The concept `StreamLinesTraits_2` describes the set of requirements to be 
fulfilled by any class used to instantiate the template parameter of the class 
`Regular_grid_2<StreamLinesTraits_2>`. 
This concept provides the types handled by the 
`Stream_lines_2<VectorField_2, Integrator_2>` class. 

\hasModel The kernels of \cgal are models for this traits class. 

*/
class StreamLinesTraits_2 {
public:

/// \name Types 
/// @{

/*! 
The scalar type. 
*/ 
typedef Hidden_type FT; 

/*! 
The point type. 
*/ 
typedef Hidden_type Point_2; 

/*! 
The vector type. 
*/ 
typedef Hidden_type Vector_2; 

/// @}

}; /* end StreamLinesTraits_2 */

