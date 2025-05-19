namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2Ref

The `Greene_convex_decomposition_2` class implements the approximation algorithm of
Greene for the decomposition of an input polygon into convex
sub-polygons \cgalCite{g-dpcp-83}. This algorithm takes \cgalBigO{n \log n}
time and \cgalBigO{n} space, where \f$ n\f$ is the size of the input polygon,
and outputs a decomposition whose size is guaranteed to be no more
than four times the size of the optimal decomposition.

\tparam Kernel must be a geometric kernel that can be used for the polygon.
\tparam Container must be a container that can be used for the polygon.
It is by default `std::vector<typename Kernel::Point_2>`.

\cgalModels{PolygonConvexDecomposition_2}

\sa `CGAL::greene_approx_convex_partition_2()`

*/
template< typename Kernel, typename Container >
class Greene_convex_decomposition_2 {
public:

  typedef CGAL::Polygon_2<Kernel,Container> Polygon_2;

}; /* end Greene_convex_decomposition_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2Ref

The `Hertel_Mehlhorn_convex_decomposition_2` class implements the approximation algorithm of Hertel
and Mehlhorn for decomposing a polygon into convex
sub-polygons \cgalCite{hm-ftsp-83}. This algorithm constructs a
triangulation of the input polygon and proceeds by removing
unnecessary triangulation edges. Given the triangulation, the
algorithm requires \cgalBigO{n} time and space to construct a convex
decomposition (where \f$ n\f$ is the size of the input polygon), whose
size is guaranteed to be no more than four times the size of the
optimal decomposition.

\tparam Kernel must be a geometric kernel that can be used for the polygon.
\tparam Container must be a container that can be used for the polygon.
It is by default `std::vector<typename Kernel::Point_2>`.

\cgalModels{PolygonConvexDecomposition_2}

\sa `CGAL::approx_convex_partition_2()`

*/
template< typename Kernel, typename Container >
class Hertel_Mehlhorn_convex_decomposition_2 {
public:

  typedef CGAL::Polygon_2<Kernel,Container> Polygon_2;

}; /* end Hertel_Mehlhorn_convex_decomposition_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2Ref

The `Optimal_convex_decomposition_2` class provides an implementation of Greene's
dynamic programming algorithm for optimal decomposition of a
polygon into convex sub-polygons \cgalCite{g-dpcp-83}. Note that
this algorithm requires \cgalBigO{n^4} time and \cgalBigO{n^3} space in
the worst case, where \f$ n\f$ is the size of the input polygon.


\tparam Kernel must be a geometric kernel that can be used for the polygon.
\tparam Container must be a container that can be used for the polygon.
It is by default `std::vector<typename Kernel::Point_2>`.

\cgalModels{PolygonConvexDecomposition_2}

\sa `CGAL::optimal_convex_partition_2()`

*/
template< typename Kernel, typename Container >
class Optimal_convex_decomposition_2 {
public:

  typedef CGAL::Polygon_2<Kernel,Container> Polygon_2;

}; /* end Optimal_convex_decomposition_2 */
} /* end namespace CGAL */
