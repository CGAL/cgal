
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Contained_in_affine_hull_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
determines whether \f$ p\f$ is contained in the 
affine hull of the points in `A = tuple [first,last)`. 
\pre The objects are of the same dimension. 

\tparam ForwardIterator has `Kernel_d::Point_d` as value type.
*/ 
template <class ForwardIterator> Bounded_side 
operator()( ForwardIterator first, ForwardIterator last, const 
Kernel_d::Point_d& p); 

/// @}

}; /* end Kernel_d::Contained_in_affine_hull_d */

