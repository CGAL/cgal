
/*!
\ingroup PkgCircularKernel2Concepts
\cgalconcept

A model `fo` of this type must provide: 

*/

class AlgebraicKernelForCircles::Solve {
public:

/// \name Operations
/// @{

/*! 
Copies in the output iterator the common roots of `p1` and `p2`, 
with their multiplicity, as objects of type 
`std::pair< AlgebraicKernelForCircles::Root_for_circles_2_2, int>`. 
*/ 
template < class OutputIterator > 
OutputIterator 
operator()(const AlgebraicKernelForCircles::Polynomial_1_2 &p1, 
const AlgebraicKernelForCircles::Polynomial_1_2 &p2, 
OutputIterator res); 

/*! 
Same as previous. 
*/ 
template < class OutputIterator > 
OutputIterator 
operator()(const AlgebraicKernelForCircles::Polynomial_1_2 &p1, 
const AlgebraicKernelForCircles::Polynomial_for_circles_2_2 &p2, 
OutputIterator res); 

/*! 
Same as previous. 
*/ 
template < class OutputIterator > 
OutputIterator 
operator()(const AlgebraicKernelForCircles::Polynomial_for_circles_2_2 &p1, 
const AlgebraicKernelForCircles::Polynomial_1_2 &p2, 
OutputIterator res); 

/*! 
Same as previous. 
*/ 
template < class OutputIterator > 
OutputIterator 
operator()(const AlgebraicKernelForCircles::Polynomial_for_circles_2_2 &p1, 
const AlgebraicKernelForCircles::Polynomial_for_circles_2_2 &p2, 
OutputIterator res); 

/// @}

}; /* end AlgebraicKernelForCircles::Solve */

