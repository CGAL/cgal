
/*!
\ingroup PkgCircularKernel2AlgebraicConcepts
\cgalConcept

\sa `CircularKernel::ConstructCircle_2`
\sa `CircularKernel::GetEquation`

*/

class AlgebraicKernelForCircles::ConstructPolynomialForCircles_2_2 {
public:

/// \name Operations
/// A model of this concept must provide:
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

