
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Linearly_independent_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*! 
decides 
whether the vectors in `A = tuple [first,last)` are linearly 
independent.

\pre The objects in `A` are of the same dimension. 
\cgalRequires The value type of `ForwardIterator` is `Kernel_d::Vector_d`. 
*/ 
template <class ForwardIterator> bool 
operator()(ForwardIterator first, ForwardIterator last); 

/// @}

}; /* end Kernel_d::Linearly_independent_d */

