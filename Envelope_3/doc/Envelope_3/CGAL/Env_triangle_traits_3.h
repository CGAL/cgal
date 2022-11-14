namespace CGAL {

/*!
\ingroup PkgEnvelope3Ref

The traits class `Env_triangle_traits_3` models the `EnvelopeTraits_3`
concept, and is used for the construction of lower and upper envelopes
of triangles in the space. It is parameterized by a \cgal-kernel,
which is parameterized in turn by a number type. The number type should
support exact rational arithmetic, to avoid numerical errors and
robustness problems. In particular, the number type should support
the arithmetic operations \f$ +\f$, \f$ -\f$, \f$ \times\f$, and \f$ \div\f$ without loss
of precision. For optimal performance, we recommend instantiating the
traits class with the predefined
`Exact_predicates_exact_constructions_kernel` provided by \cgal.
Using this kernel guarantees exactness and robustness, while it incurs
only a minor overhead (in comparison to working with a fast, inexact number
type) for most inputs.

Note that when we project the boundary of a triangle, or the intersection
of two triangles, onto the \f$ xy\f$-plane, we obtain line segments. Indeed,
`Env_triangle_traits_3` inherits from the `Arr_segment_traits_2<Kernel>` traits
class, and extends it by adding operations on 3D objects, namely spacial
triangles. Note that the traits class does <I>not</I> define
`Kernel::Triangle_3` as its surface (and \f$ xy\f$-monotone surface) type,
as one may expect. This is because the traits class needs to store extra
data with the triangles in order to efficiently operate on them.
Nevertheless, the nested `Xy_monotone_surface_3` and `Surface_3`
types are however constructible from a `Kernel::Triangle_3`
instance and are also convertible to a `Kernel::Triangle_3` object.
Both types, `Xy_monotone_surface_3` and `Surface_3`, refer to the
same class, as <I>every</I> triangle is (weakly) \f$ xy\f$-monotone).

\cgalModels `EnvelopeTraits_3`

*/
template< typename Kernel >
class Env_triangle_traits_3 : public Arr_segment_traits_2<Kernel> {
public:

}; /* end Env_triangle_traits_3 */
} /* end namespace CGAL */
