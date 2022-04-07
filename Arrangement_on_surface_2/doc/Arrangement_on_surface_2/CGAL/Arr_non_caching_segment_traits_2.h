
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2TraitsClasses

The traits class `Arr_non_caching_segment_traits_2` is a model of the `ArrangementTraits_2`
concept that allows the construction and maintenance of arrangements of
line segments. It is parameterized with a \cgal-Kernel type, and it
is derived from it. This traits class is a thin layer above the
parameterized kernel. It inherits the `Point_2` from the kernel and its
`X_monotone_curve_2` and `Curve_2` types are both defined as
`Kernel::Segment_2`. Most traits-class functor are inherited from the
kernel functor, and the traits class only supplies the necessary functors
that are not provided by the kernel. The kernel is parameterized with a
number type, which should support exact rational arithmetic in order to
avoid robustness problems, although other number types could be used at the
user's own risk.

The traits-class implementation is very simple, yet may lead to
a cascaded representation of intersection points with exponentially long
bit-lengths, especially if the kernel is parameterized with a number type
that does not perform normalization (e.g. `Quotient<MP_Float>`).
The `Arr_segment_traits_2` traits class avoids this cascading
problem, and should be the default choice for implementing arrangements of
line segments. It is recommended to use `Arr_non_caching_segment_traits_2` only for very sparse
arrangements of huge sets of input segments.

While `Arr_non_caching_segment_traits_2` models the concept
`ArrangementDirectionalXMonotoneTraits_2`, the implementation of
the `Are_mergeable_2` operation does not enforce the input curves
to have the same direction as a precondition. Moreover, `Arr_non_caching_segment_traits_2`
supports the merging of curves of opposite directions.

\cgalModels `ArrangementTraits_2`
\cgalModels `ArrangementLandmarkTraits_2`
\cgalModels `ArrangementDirectionalXMonotoneTraits_2`

\sa `Arr_segment_traits_2<Kernel>`

*/
template< typename Kernel >
class Arr_non_caching_segment_traits_2
  : public Arr_non_caching_segment_basic_traits_2<Kernel>
 {
public:

}; /* end Arr_non_caching_segment_traits_2 */
} /* end namespace CGAL */
