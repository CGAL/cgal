
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

\sa `CircularKernel::ConstructLine_2`
\sa `CircularKernel::GetEquation`

*/

class AlgebraicKernelForCircles::ConstructPolynomial_1_2 {
public:

/// \name Operations
/// @{

/*! 
Constructs polynomial `ax+by+c`. 
*/ 
AlgebraicKernelForCircles::Polynomial_1_2 
operator()(const AlgebraicKernelForCircles::RT &a, 
const AlgebraicKernelForCircles::RT &b, 
const AlgebraicKernelForCircles::RT &c); 

/// @}

}; /* end AlgebraicKernelForCircles::ConstructPolynomial_1_2 */

