
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

\refines `Kernel::ConstructLine_2`

\sa `CircularKernel::GetEquation`

*/

class CircularKernel::ConstructLine_2 {
public:

/// \name Operations
///  A model `fo` of this type must provide:
/// @{

/*! 
Constructs the supporting line of a line segment. 
*/ 
CircularKernel::Line_2 operator() 
(CircularKernel::LineArc_2); 

/*! 
Constructs a line from an equation. 
*/ 
CircularKernel::Line_2 operator() 
(CircularKernel::Polynomial_1_2); 

/// @}

}; /* end CircularKernel::ConstructLine_2 */

