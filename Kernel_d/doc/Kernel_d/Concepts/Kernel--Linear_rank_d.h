
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

*/

class Kernel_d::Linear_rank_d {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
computes
the linear rank of the vectors in `A = tuple [first,last)`.
\pre \f$ A\f$ contains vectors of the same dimension \f$ d\f$.

\tparam ForwardIterator has `Kernel_d::Vector_d` as value type.
*/
template <class ForwardIterator> int
operator()(ForwardIterator first, ForwardIterator last);

/// @}

}; /* end Kernel_d::Linear_rank_d */

