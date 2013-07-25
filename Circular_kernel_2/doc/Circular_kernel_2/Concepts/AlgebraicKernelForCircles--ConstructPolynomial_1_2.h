
/*!
\ingroup PkgCircularKernel2AlgebraicConcepts
\cgalConcept

\sa `CircularKernel::ConstructLine_2`
\sa `CircularKernel::GetEquation`

*/

class AlgebraicKernelForCircles::ConstructPolynomial_1_2 {
public:

/// \name Operations
/// A model of this concept must provide: 
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

