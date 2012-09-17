
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

\sa `CircularKernel::ConstructCircle_2`
\sa `CircularKernel::GetEquation`

*/

class AlgebraicKernelForCircles::ConstructPolynomialForCircles_2_2 {
public:

/// \name See Also 
/// @{

/*! 
Constructs polynomial `(x-a)^2 + (y-b)^2 - rsq`. 
*/ 
AlgebraicKernelForCircles::PolynomialForCircles_2_2 
operator()(const AlgebraicKernelForCircles::FT a, 
const AlgebraicKernelForCircles::FT b, 
const AlgebraicKernelForCircles::FT rsq); 

/// @}

}; /* end AlgebraicKernelForCircles::ConstructPolynomialForCircles_2_2 */

