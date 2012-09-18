
/*!
\ingroup PkgKernelDKernelConcept
\cgalconcept

A model for this must provide: 

*/

class Kernel_d::Linear_base_d {
public:

/// \name See Also 
/// @{

/*! 
computes a basis of the linear space 
spanned by the vectors in `A = tuple [first,last)` and returns 
it via an iterator range starting in `result`. The returned 
iterator marks the end of the output.

\pre \f$ A\f$ contains vectors of the same dimension \f$ d\f$. 
\requires The value type of `ForwardIterator` and `OutputIterator` is `Kernel_d::Vector_d`. 
*/ 
template <class 
ForwardIterator, class OutputIterator> int 
operator()(ForwardIterator first, ForwardIterator last, 
OutputIterator result); 

/// @}

}; /* end Kernel_d::Linear_base_d */

