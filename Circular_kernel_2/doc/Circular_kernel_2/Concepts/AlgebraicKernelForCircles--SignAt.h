
/*!
\ingroup PkgCircularKernel2AlgebraicConcepts
\cgalConcept

*/

class AlgebraicKernelForCircles::SignAt {
public:

/// \name Operations
/// A model of this concept must provide: 
/// @{

/*!
Computes the sign of polynomial `p` evaluated at a root `r`. 
*/ 
template < class OutputIterator > 
CGAL::Sign 
operator()(const AlgebraicKernelForCircles::Polynomial_1_2 & p, 
const AlgebraicKernelForCircles::Root_for_circles_2_2 & r); 

/*!
Same as previous. 
*/ 
template < class OutputIterator > 
CGAL::Sign 
operator()(const AlgebraicKernelForCircles::Polynomial_for_circles_2_2 & p, 
const AlgebraicKernelForCircles::Root_for_circles_2_2 & r); 

/// @}

}; /* end AlgebraicKernelForCircles::SignAt */

