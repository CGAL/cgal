namespace CGAL {

/*!
\ingroup PkgEnvelope3Ref

The traits class `Env_plane_traits_3` models the `EnvelopeTraits_3` concept,
and is used for the construction of lower and upper envelopes of planes
and half planes in the space. It is parameterized by a \cgal-kernel model,
which is parameterized in turn by a number type. The number type should
support exact rational arithmetic, to avoid numerical errors and
robustness problems. In particular, the number type should support the
arithmetic operations \f$ +\f$, \f$ -\f$, \f$ \times\f$, and \f$ \div\f$ without loss of
precision. For optimal performance, we recommend instantiating the traits
class with the predefined
`Exact_predicates_exact_constructions_kernel` provided by \cgal.
Using this kernel guarantees exactness and robustness, while it incurs
only a minor overhead (in comparison to working with a fast, inexact number
type) for most inputs.

Note that an entire plane has no boundaries, and the projection of a
half-plane is an (unbounded) line. Naturally, rays and segments may occur as
a result of overlaying projections of several half planes. Indeed,
`Env_plane_traits_3` inherits from the `Arr_linear_traits_2<Kernel>` traits
class, and extends it by adding operations on planes and half planes.
The nested `Xy_monotone_surface_3` and `Surface_3` types refer
to the same type. They are constructible from a `Kernel::Plane_3`
in case of an entire plane, or from `Kernel::Plane_3` and
`Kernel::Line_2` in case of a half-plane. The line orientation
determines which half is considered.

\cgalModels `EnvelopeTraits_3`

*/
template< typename Kernel >
class Env_plane_traits_3 : public Arr_linear_traits_2<Kernel> {
public:

}; /* end Env_plane_traits_3 */
} /* end namespace CGAL */
