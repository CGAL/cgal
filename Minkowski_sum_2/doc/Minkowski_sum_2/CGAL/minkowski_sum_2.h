namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons.
This method defaults to the reduced convolution method, see below.
Note that as the input polygons may not be convex, their Minkowski
sum may not be a simple polygon. The result is therefore represented
as a polygon with holes.

\pre Both `P` and `Q` are simple, counterclockwise-oriented polygons.

\sa `CGAL::minkowski_sum_reduced_convolution_2()`
\sa `CGAL::minkowski_sum_full_convolution_2()`
 */
template<class Kernel, class Container>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_2 (const Polygon_2<Kernel,Container>& P,
const Polygon_2<Kernel,Container>& Q);

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons. The
function computes the reduced convolution \cgalCite{cgal:bl-frmsurc-11} of
the two polygons and extracts those loops of the convolution which are part of
the Minkowsi sum. This method works very efficiently, regardless of whether `P`
and `Q` are convex or non-convex. It is usually faster than the full
convolution method, except in degenerate cases where the output polygon has
many holes.
Note that as the input polygons may not be convex, their Minkowski
sum may not be a simple polygon. The result is therefore represented
as a polygon with holes.

\pre Both `P` and `Q` are simple, counterclockwise-oriented polygons.
*/

template<class Kernel, class Container>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_reduced_convolution_2 (const Polygon_2<Kernel,Container>& P,
const Polygon_2<Kernel,Container>& Q);

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons.
The function computes the convolution cycles of the two polygons and
extract the regions having positive winding number with respect to these
cycles. This method work very efficiently, regardless of whether `P`
and `Q` are convex or non-convex.
Note that as the input polygons may not be convex, their Minkowski
sum may not be a simple polygon. The result is therefore represented
as a polygon with holes.

\pre Both `P` and `Q` are simple polygons.
*/
template<class Kernel, class Container>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_full_convolution_2 (const Polygon_2<Kernel,Container>& P,
const Polygon_2<Kernel,Container>& Q);

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons.
If the input polygons `P` and `Q` are not convex, the function
decomposes them into convex sub-polygons \f$ P_1, \ldots, P_k\f$ and
\f$ Q_1, \ldots, Q_{\ell}\f$ and computes the union of pairwise sub-sums
(namely \f$ \bigcup_{i,j}{(P_i \oplus Q_j)}\f$).
The decomposition is performed using the given decomposition strategy
`decomp`, which must be an instance of a class that models the
concept `PolygonConvexDecomposition`.
Note that as the input polygons may not be convex, their Minkowski
sum may not be a simple polygon. The result is therefore represented
as a polygon with holes.
\pre Both `P` and `Q` are simple polygons.
*/
template<class Kernel, class Container,
class DecompositionStrategy>
Polygon_with_holes_2<Kernel,Container>
minkowski_sum_2 (const Polygon_2<Kernel,Container>& P,
const Polygon_2<Kernel,Container>& Q,
const DecompositionStrategy& decomp);

} /* namespace CGAL */
