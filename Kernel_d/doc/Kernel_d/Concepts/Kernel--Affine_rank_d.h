
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Affine_rank_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*! 
computes 
the affine rank of the points in `A = tuple [first,last)`. 
\pre The objects are of the same dimension. 

\cgalRequires The value type of `ForwardIterator` is `Kernel_d::Point_d`. 
*/ 
template <class ForwardIterator> int 
operator()(ForwardIterator first, ForwardIterator last); 

/// @}

}; /* end Kernel_d::Affine_rank_d */

