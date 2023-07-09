
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2TraitsClasses

The traits class `Arr_non_caching_segment_basic_traits_2` is a model of the `ArrangementTraits_2`
concept that allow the construction and maintenance of arrangements of
sets of pairwise interior-disjoint line segments. It is templated with a
\cgal-Kernel model, and it is derived from it. This traits class is a
thin layer above the parameterized kernel. It inherits the `Point_2`
from the kernel and its `X_monotone_curve_2` type is defined as
`Kernel::Segment_2`. Most traits-class functor are inherited from the
kernel functor, and the traits class only supplies the necessary functors
that are not provided by the kernel. The kernel is parameterized with a
number type, which should support the arithmetic operations \f$ +\f$, \f$ -\f$ and
\f$ \times\f$ in an exact manner in order to avoid robustness problems.
Using `Cartesian<MP_Float>` or `Cartesian<Gmpz>` are possible
instantiations for the kernel. Using other (inexact) number types
(for example, instantiating the template with
`Simple_cartesian<double>`) is also possible, at the user's own
risk.

\cgalModels{ArrangementLandmarkTraits_2}

*/
template< typename Kernel >
class Arr_non_caching_segment_basic_traits_2 {
public:

}; /* end Arr_non_caching_segment_basic_traits_2 */
} /* end namespace CGAL */
