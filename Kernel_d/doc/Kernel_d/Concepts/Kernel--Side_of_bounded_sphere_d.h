
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Side_of_bounded_sphere_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
returns the relative position of point 
`p` to the sphere defined by `A = tuple [first,last)`. The 
order of the points of \f$ A\f$ does not matter.

\pre `orientation(first,last)` is not `ZERO`. 
\tparam ForwardIterator has `Kernel_d::Point_d` as value type.
*/ 
template <class ForwardIterator> Bounded_side 
operator()( ForwardIterator first, ForwardIterator last, const 
Kernel_d::Point_d& p); 

/// @}

}; /* end Kernel_d::Side_of_bounded_sphere_d */

