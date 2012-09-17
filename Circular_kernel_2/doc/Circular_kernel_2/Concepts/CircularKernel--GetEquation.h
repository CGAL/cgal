
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

\sa `CircularKernel::ConstructLine_2`
\sa `CircularKernel::ConstructCircle_2`

*/

class CircularKernel::GetEquation {
public:

/// \name Operations 
/// @{

/*! 
Returns the equation of the line. 
*/ 
CircularKernel::Polynomial_1_2 
operator()(const CircularKernel::Line_2 & c); 

/*! 
Returns the equation of the circle. 
*/ 
CircularKernel::Polynomial_for_circles_2_2 
operator()(const CircularKernel::Circle_2 & c); 

/// @}

}; /* end CircularKernel::GetEquation */

