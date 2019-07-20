namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2Ref

Provides a guaranteed approximation of the inset, or inner offset, of
the given polygon `P` by a given radius `r`. Namely, the
function computes the set of points inside the polygon whose distance
from \f$ P\f$'s boundary is at least \f$ r\f$:
\f$ \{ p \in P \;|\; {\rm dist}(p, \partial P) \geq r \}\f$,
with the approximation error bounded by `eps`.
Note that as the input polygon may not be convex, its inset may comprise
several disconnected components. The result is therefore represented as a
sequence of generalized polygons, whose edges are either line segments or
circular arcs.
The output sequence is returned via the output iterator `oi`, whose
value-type must be `Gps_circle_segment_traits_2::Polygon_2`.
\pre `P` is a simple polygon.
*/
template<class Kernel, class Container, class OutputIterator>
OutputIterator
approximated_inset_2 (const Polygon_2<Kernel, Container>& pgn,
const typename Kernel::FT& r,
const double& eps,
OutputIterator oi);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2Ref

Provides a guaranteed approximation of the offset of the given polygon
`P` by a given radius `r` - namely, the function computes the
Minkowski sum \f$ P \oplus B_r\f$, where \f$ B_r\f$ is a disc of radius
`r` centered at the origin.
The function actually outputs a set \f$ S\f$ that contains the Minkowski sum,
such that the approximation error is bounded by `eps`.
Note that as the input polygon may not be convex, its offset may not be a
simple polygon. The result is therefore represented as a polygon with
holes, whose edges are either line segments or circular arcs.
\pre `P` is a simple polygon.
*/
template<class Kernel, class Container>
typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
approximated_offset_2 (const Polygon_2<Kernel, Container>& P,
const typename Kernel::FT& r,
const double& eps);

/*!
\ingroup PkgMinkowskiSum2Ref

Provides a guaranteed approximation of offset the given polygon with holes
`pwh` by a given radius `r`, such that the approximation error is bounded
by `eps`. It does so by offsetting outer boundary of `pwh` and insetting
its holes.
The result is represented as a generalized polygon with holes, such that the edges
of the polygon correspond to line segment and circular arcs.
\pre `pwh` is <I>not</I> unbounded (it has a valid outer boundary).
*/
template<class Kernel, class Container>
typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
approximated_offset_2 (const Polygon_with_holes_2<Kernel, Container>& wh,
const typename Kernel::FT& r,
const double& eps);

/*!
\ingroup PkgMinkowskiSum2Ref

Provides a guaranteed approximation of the offset of the given polygon
`P` by a radius `r`, as described above.
If the input polygon `P` is not convex, the function
decomposes it into convex sub-polygons \f$ P_1, \ldots, P_k\f$ and computes
the union of the sub-offsets (namely \f$ \bigcup_{i}{(P_i \oplus B_r)}\f$).
The decomposition is performed using the given decomposition strategy
`decomp`, which must be an instance of a class that models the
concept `PolygonConvexDecomposition`.
\pre `P` is a simple polygon.
*/
template<class Kernel, class Container,
class DecompositionStrategy>
typename Gps_circle_segment_traits_2<Kernel>::Polygon_with_holes_2
approximated_offset_2 (const Polygon_2<Kernel, Container>& P,
const typename Kernel::FT& r,
const double& eps,
const DecompositionStrategy& decomp);

} /* namespace CGAL */
