
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Contained_in_simplex_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*! 
determines whether \f$ p\f$ is contained in the 
simplex of the points in `A = tuple [first,last)`.

\pre The objects in \f$ A\f$ are of the same dimension and affinely independent. 
\cgalRequires The value type of `ForwardIterator` is `Kernel_d::Point_d`. 
*/ 
template <class ForwardIterator> Bounded_side 
operator()( ForwardIterator first, ForwardIterator last, const 
Kernel_d::Point_d& p); 

/// @}

}; /* end Kernel_d::Contained_in_simplex_d */

