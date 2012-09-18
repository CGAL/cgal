
/*!
\ingroup PkgKernelDKernelConcept
\cgalconcept

A model for this must provide: 

*/

class Kernel_d::Affinely_independent_d {
public:

/// \name Has Models 
/// @{

/*! 
returns 
true iff the points in `A = tuple [first,last)` are affinely 
independent.

\pre The objects are of the same dimension. 

\requires The value type of `ForwardIterator` is `Kernel_d::Point_d`. 
*/ 
template <class ForwardIterator> bool 
operator()(ForwardIterator first, ForwardIterator last); 

/// @}

}; /* end Kernel_d::Affinely_independent_d */

