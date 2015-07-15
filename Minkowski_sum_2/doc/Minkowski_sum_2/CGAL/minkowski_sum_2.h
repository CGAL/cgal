namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of two given polygons
(which may have holes). `PolygonType1` and `PolygonType2` can be any combination of:

- `Polygon_2`
- `Polygon_with_holes_2`

This method defaults to the reduced convolution method, see below.

\sa `CGAL::minkowski_sum_by_reduced_convolution_2()`
\sa `CGAL::minkowski_sum_by_full_convolution_2()`
 */
template <typename Kernel, typename Container>
Polygon_with_holes_2<Kernel, Container>
minkowski_sum_2(const PolygonType1<Kernel, Container>& P,
                const PolygonType2<Kernel, Container>& Q);

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons. The
function computes the reduced convolution \cgalCite{cgal:bl-frmsurc-11} of
the two polygons and extracts those loops of the convolution that are part of
the Minkowsi sum. This method works very efficiently, regardless of whether `P`
and `Q` are convex or non-convex. In general, it is more efficient than the
full convolution method, except for degenerate cases where the output polygon
has many holes.

`PolygonType1` and `PolygonType2` can be any combination of:

- `Polygon_2`
- `Polygon_with_holes_2`
*/

template <typename Kernel, typename Container>
Polygon_with_holes_2<Kernel, Container>
minkowski_sum_by_reduced_convolution_2(const PolygonType1<Kernel, Container>& P,
                                       const PolygonType2<Kernel, Container>& Q);

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons.
The function computes the (full) convolution cycles of the two polygons and
extract the regions having positive winding number with respect to these
cycles. This method work very efficiently, regardless of whether `P`
and `Q` are convex or non-convex.

`PolygonType1` and `PolygonType2` can be any combination of:

- `Polygon_2`
- `Polygon_with_holes_2`
*/
template <typename Kernel, typename Container>
Polygon_with_holes_2<Kernel, Container>
minkowski_sum_by_full_convolution_2(const Polygon_2<Kernel, Container>& P,
                                    const Polygon_2<Kernel, Container>& Q);

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons.
If the input polygons `P` and `Q` are not convex, the function
decomposes them into convex sub-polygons \f$ P_1, \ldots, P_k\f$ and
\f$ Q_1, \ldots, Q_{\ell}\f$ and computes the union of pairwise sub-sums
(namely \f$ \bigcup_{i,j}{(P_i \oplus Q_j)}\f$).
The decomposition is performed using the given decomposition method
`decomp`, which must be an instance of a class template that models
the concept `PolygonConvexDecomposition_2`.
*/
template <typename Kernel, typename Container, typename PolygonConvexDecomposition_2>
Polygon_with_holes_2<Kernel, Container>
minkowski_sum_2(const PolygonType1<Kernel, Container>& P,
                const PolygonType2<Kernel, Container>& Q,
                const PolygonConvexDecomposition_2& decomp,
                const Gps_segment_traits_2& traits = Gps_segment_traits_2<Kernel,Container, Arr_segment_traits>());

} /* namespace CGAL */
