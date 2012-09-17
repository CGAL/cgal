
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

*/

class AlgebraicKernelForCircles::SignAt {
public:

/// \name Operations
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

