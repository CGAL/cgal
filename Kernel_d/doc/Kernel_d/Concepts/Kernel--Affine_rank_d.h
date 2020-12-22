
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

\tparam ForwardIterator is a model of `ForwardIterator` with `Kernel_d::Point_d` as value type.
*/
template <class ForwardIterator> int
operator()(ForwardIterator first, ForwardIterator last);

/// @}

}; /* end Kernel_d::Affine_rank_d */

