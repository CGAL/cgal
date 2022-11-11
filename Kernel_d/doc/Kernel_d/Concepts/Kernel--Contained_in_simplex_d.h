
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Contained_in_simplex_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
determines whether \f$ p\f$ is contained in the
simplex of the points in `A = tuple [first,last)`.

\pre The objects in \f$ A\f$ are of the same dimension and affinely independent.
\tparam ForwardIterator has `Kernel_d::Point_d` as value type.
*/
template <class ForwardIterator> Bounded_side
operator()( ForwardIterator first, ForwardIterator last, const
Kernel_d::Point_d& p);

/// @}

}; /* end Kernel_d::Contained_in_simplex_d */

