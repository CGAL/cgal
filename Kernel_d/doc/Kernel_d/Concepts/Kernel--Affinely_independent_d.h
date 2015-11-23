
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept
*/

class Kernel_d::Affinely_independent_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns 
true iff the points in `A = tuple [first,last)` are affinely 
independent.

\pre The objects are of the same dimension. 

\tparam ForwardIterator has `Kernel_d::Point_d` as value type.

*/ 
template <class ForwardIterator> bool 
operator()(ForwardIterator first, ForwardIterator last); 

/// @}

}; /* end Kernel_d::Affinely_independent_d */

