
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
\tparam ForwardIterator has `Kernel_d::Vector_d` as value type.
*/
template <class ForwardIterator> Bounded_side
operator()( ForwardIterator first, ForwardIterator last, const
Kernel_d::Vector_d& v);

/// @}

}; /* end Kernel_d::Contained_in_linear_hull_d */

