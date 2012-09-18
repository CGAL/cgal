
/*!
\ingroup PkgKernelDKernelConcept
\cgalconcept

A model for this must provide: 

*/

class Kernel_d::Affine_rank_d {
public:

/// \name Has Models 
/// @{

/*! 
computes 
the affine rank of the points in `A = tuple [first,last)`. 
\pre The objects are of the same dimension. 

\requires The value type of `ForwardIterator` is `Kernel_d::Point_d`. 
*/ 
template <class ForwardIterator> int 
operator()(ForwardIterator first, ForwardIterator last); 

/// @}

}; /* end Kernel_d::Affine_rank_d */

