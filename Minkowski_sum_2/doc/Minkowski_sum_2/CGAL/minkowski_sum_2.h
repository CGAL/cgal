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
The function decomposes the summand `P` into convex sub-polygons
\f$ P_1, \ldots, P_k\f$ using the given decomposition method `decomp_P`.
If the summand `P` is of type `Polygon_2`, then `decomp_P` must be an
instance of a class template that models the concept
`PolygonConvexDecomposition_2`. If `P` is of type `Polygon_with_holes_2`,
then `decomp_P` must be an instance of a class template that models the
concept `PolygonWithHolesConvexDecomposition_2`.
Similarly, the function decomposes the summand `Q` into convex sub-polygons
\f$ Q_1, \ldots, Q_k\f$ using the given decomposition method `decomp_Q`.
If the summand `Q` is of type `Polygon_2`, then `decomp_Q` must be an
instance of a class template that models the concept
`PolygonConvexDecomposition_2`. If `Q` is of type `Polygon_with_holes_2`,
then `decomp_Q` must be an instance of a class template that models the
concept `PolygonWithHolesConvexDecomposition_2`.
Then, the function computes the union of pairwise sub-sums (namely
\f$ \bigcup_{i,j}{(P_i \oplus Q_j)}\f$).
*/
template <typename Kernel, typename Container,
          typename PolygonConvexDecompositionP_2_,
          typename PolygonConvexDecompositionQ_2_>
Polygon_with_holes_2<Kernel, Container>
minkowski_sum_2(const PolygonType1<Kernel, Container>& P,
                const PolygonType2<Kernel, Container>& Q,
                const PolygonConvexDecompositionP_2_& decomp_P,
                const PolygonConvexDecompositionQ_2_& decomp_Q,
                const Gps_segment_traits_2& traits = Gps_segment_traits_2<Kernel,Container, Arr_segment_traits>());

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons.
The function decomposes the summands `P` and `Q` into convex sub-polygons
\f$ P_1, \ldots, P_k\f$ and \f$ Q_1, \ldots, Q_{\ell}\f$, respectively.
Then, it computes the union of pairwise sub-sums (namely
\f$ \bigcup_{i,j}{(P_i \oplus Q_j)}\f$). The decomposition is performed
using the given decomposition method `decomp`. If both summands
\f$ P \oplus Q\f$ are of type `Polygon_2`, then `decomp` must be instance
of a class template that models the concept `PolygonConvexDecomposition_2`.
If both summands are of type `Polygon_with_holes_2`, then `decomp` must be
a model of the concept `PolygonWithHolesConvexDecomposition_2`.
*/
template <typename Kernel, typename Container,
          typename PolygonConvexDecomposition_2_>
Polygon_with_holes_2<Kernel, Container>
minkowski_sum_2(const PolygonType1<Kernel, Container>& P,
                const PolygonType2<Kernel, Container>& Q,
                const PolygonConvexDecomposition_2_& decomp,
                const Gps_segment_traits_2& traits = Gps_segment_traits_2<Kernel,Container, Arr_segment_traits>());

/*!
\ingroup PkgMinkowskiSum2

Computes the Minkowski sum \f$ P \oplus Q\f$ of the two given polygons
using the decomposition strategy. It decomposes the summands `P` and `Q`
into convex sub-polygons \f$ P_1, \ldots, P_k\f$ and
\f$ Q_1, \ldots, Q_{\ell}\f$, respectively. Then, it computes the union
of pairwise sub-sums (namely \f$ \bigcup_{i,j}{(P_i \oplus Q_j)}\f$).

If the summand `P` is of type `Polygon_with_holes_2`, then the function
first applies the hole filteration on `P`. If the summand `P` remains a
polygon with holes, then the function decomposes the summand `P` using
the given decomposition method `with_holes_decomp`. If, however, `P`
turns into a polygon without holes, then the function decomposes the
summand `P` using the given decomposition method `no_holes_decomp`,
unless the result is a convex polygon, in which case the nop strategy is
applied; namely, an instance of the class template
`Polygon_nop_decomposition_2` is used. If `P` is a polygon without holes
to start with, then only convexity is checked (checking whether the
result is convex inccurs a small overhead though). Then depending on
the result either `no_holes_decomp` or the nop strategy is applied.
Similarly, if the summand `Q` is of type `Polygon_with_holes_2`, then
the function first applies the hole filteration on `Q`. If the summand
`Q` remains a polygon with holes, then the function decomposes the
summand `Q` using the given decomposition method `with_holes_decomp`.
If, however, `Q` turns into a polygon without holes, then the function
decomposes the summand `Q` using the given decomposition method
`no_holes_decomp`, unless the result is a convex polygon, in which case
the nop strategy is applied. If `Q` is a polygon without holes to start
with, then only convexity is checked and the decomposition strategy is
chosen accordingly.

\tparam PolygonNoHolesConvexDecomposition_2_ a model of the concept `PolygonConvexDecomposition_2`.
\tparam PolygonWithHolesConvexDecomposition_2_ a model of the concept `PolygonWithHolesConvexDecomposition_2`.
*/
  template <typename Kernel, typename Container,
            typename PolygonNoHolesConvexDecomposition_2_,
            typename PolygonWithHolesConvexDecomposition_2_>
Polygon_with_holes_2<Kernel, Container>
minkowski_sum_by_decomposition_2(const PolygonType1<Kernel, Container>& P,
                const PolygonType2<Kernel, Container>& Q,
                const PolygonNoHolesConvexDecomposition_2_& no_holes_decomp,
                const PolygonWithHolesConvexDecomposition_2_& with_holes_decomp,
                const Gps_segment_traits_2& traits = Gps_segment_traits_2<Kernel,Container, Arr_segment_traits>());

} /* namespace CGAL */
