namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2

Computes the inset, or inner offset, of the given polygon `P` by a
given radius `r` - namely, the function computes the set of points
inside the polygon whose distance from \f$ P\f$'s boundary is at least \f$ r\f$:
\f$ \{ p \in P \;|\; {\rm dist}(p, \partial P) \geq r \}\f$.
Note that as the input polygon may not be convex, its inset may comprise
several disconnected components. The result is therefore represented as a
sequence of generalized polygons, such that the edges of each polygon
correspond to line segments and circular arcs, both are special types of
conic arcs, as represented by the `traits` class.
The output sequence is returned via the output iterator `oi`, whose
value-type must be `Gps_traits_2::Polygon_2`.
\pre `P` is a simple polygon.
*/
template<class ConicTraits, class Container, class OutputIterator>
OutputIterator inset_polygon_2
(const Polygon_2<typename ConicTraits::Rat_kernel,
Container>& P,
const typename ConicTraits::Rat_kernel::FT& r,
const ConicTraits& traits,
OutputIterator oi);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2

Computes the offset of the given polygon `P` by a given radius
`r` - namely, the function computes the Minkowski sum
\f$ P \oplus B_r\f$, where \f$ B_r\f$ is a disc of radius `r` centered at the
origin.
Note that as the input polygon may not be convex, its offset may not be a
simple polygon. The result is therefore represented as a generalized
polygon with holes, such that the edges of the polygon correspond to
line segments and circular arcs, both are special types of conic arcs,
as represented by the `traits` class.
\pre `P` is a simple polygon.
*/
template<class ConicTraits, class Container>
typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
offset_polygon_2
(const Polygon_2<typename ConicTraits::Rat_kernel,
Container>& P,
const typename ConicTraits::Rat_kernel::FT& r,
const ConicTraits& traits);

/*!
\ingroup PkgMinkowskiSum2

Computes the offset of the given polygon with holes `pwh` by a given
radius `r`. It does so by offsetting outer boundary of `pwh` and
insetting its holes.
The result is represented as a generalized polygon with holes, such that the
edges of the polygon correspond to line segments and circular arcs, both are
special types of conic arcs, as represented by the `traits` class.
\pre `pwh` is <I>not</I> unbounded (it has a valid outer boundary).
*/
template<class ConicTraits, class Container>
typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
offset_polygon_2
(const Polygon_with_holes_2<typename ConicTraits::Rat_kernel,
Container>& pwh,
const typename ConicTraits::Rat_kernel::FT& r,
const ConicTraits& traits);

/*!
\ingroup PkgMinkowskiSum2

Computes the exact representation of the offset of the given polygon
`P` by a radius `r`, as described above.
If `P` is not convex, the function decomposes it into convex
sub-polygons \f$ P_1, \ldots, P_k\f$ and computes the union of sub-offsets
(namely \f$ \bigcup_{i}{(P_i \oplus B_r)}\f$).
The decomposition is performed using the given decomposition strategy
`decomp`, which must be an instance of a class that models the
concept `PolygonConvexDecomposition`.
\pre `P` is a simple polygon.
*/
template<class ConicTraits, class Container,
class DecompositionStrategy>
typename Gps_traits_2<ConicTraits>::Polygon_with_holes_2
offset_polygon_2
(const Polygon_2<typename ConicTraits::Rat_kernel,
Container>& P,
const typename ConicTraits::Rat_kernel::FT& r,
const DecompositionStrategy& decomp,
const ConicTraits& traits);

} /* namespace CGAL */
