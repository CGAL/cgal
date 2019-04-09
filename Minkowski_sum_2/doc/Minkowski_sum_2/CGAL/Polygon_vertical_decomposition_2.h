namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2Ref

The `Polygon_vertical_decomposition_2` class implements a convex
decompistion of a polygon or a polygon with holes into pseudo trapezoids
utilizing the CGAL::decompose() free function of the
\ref chapterArrangement_on_surface_2 "2D Arrangements" package.

The algorithm operates in \f$ O(n \log n)\f$ time and takes
\f$ O(n)\f$ space at the worst case, where \f$ n\f$ is the
size of the input polygon.

\cgalModels `PolygonWithHolesConvexDecomposition_2`

*/
template <typename Kernel, typename Container>
class Polygon_vertical_decomposition_2 {
public:

  /// @{
  /// @}

}; /* end Polygon_vertical_decomposition_2 */
} /* end namespace CGAL */
