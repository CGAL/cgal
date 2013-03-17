
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Contained_in_linear_hull_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*! 
determines whether \f$ v\f$ is contained in the 
linear hull of the vectors in `A = tuple [first,last)`. 
\pre The objects are of the same dimension. 
\cgalRequires The value type of `ForwardIterator` is `Kernel_d::Vector_d`. 
*/ 
template <class ForwardIterator> Bounded_side 
operator()( ForwardIterator first, ForwardIterator last, const 
Kernel_d::Vector_d& v); 

/// @}

}; /* end Kernel_d::Contained_in_linear_hull_d */

