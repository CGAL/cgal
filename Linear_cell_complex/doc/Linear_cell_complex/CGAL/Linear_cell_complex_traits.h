
namespace CGAL {

/*!
\ingroup PkgLinearCellComplexClasses

This geometric traits concept is used in the
`Linear_cell_complex` class. It can take as parameter any model of the
concept `Kernel` (for example any \cgal kernel), and defines inner
types and functors corresponding to the given dimension.

\cgalModels `LinearCellComplexTraits`

\tparam d the dimension of the kernel,
\tparam K a model of the concept `Kernel` if `d==2` or
`d==3`; a model of the concept `Kernel_d` otherwise.

There is a default template arguments for `K` which is
`CGAL::Exact_predicates_inexact_constructions_kernel`
if `d` is 2 or 3, and is `CGAL::Cartesian_d<double>`
otherwise.

Note that the default argument used for `K` when
<I>d</I>>3 does not use exact predicates because operations that
use predicates are only defined in 2D and 3D.

\sa `CGAL::Linear_cell_complex<d,d2,LCCTraits,Items,Alloc>`

*/
template< typename d, typename K >
class Linear_cell_complex_traits : public K {
public:

/// \name Constants
/// @{

/*!

*/
static unsigned int ambient_dimension = d;

/// @}

}; /* end Linear_cell_complex_traits */
} /* end namespace CGAL */
