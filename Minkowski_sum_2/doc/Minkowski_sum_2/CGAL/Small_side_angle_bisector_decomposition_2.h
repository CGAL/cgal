namespace CGAL {

/*!
\ingroup PkgMinkowskiSum2Ref

The `Small_side_angle_bisector_decomposition_2` class implements a simple yet efficient heuristic for
decomposing an input polygon into convex sub-polygons. It is based
on the algorithm suggested by Flato and Halperin \cgalCite{fh-recpm-00},
but without introducing Steiner points. The algorithm operates in two
major steps. In the first step, it tries to subdivide the polygon by
connect two reflex vertices with an edge. When this is not possible any
more, it eliminates the reflex vertices one by one by connecting them
to other convex vertices, such that the new edge best approximates
the angle bisector of the reflex vertex. The algorithm operates in
\f$ O(n^2)\f$ time and takes \f$ O(n)\f$ space at the worst case, where
\f$ n\f$ is the size of the input polygon.

\cgalModels `PolygonConvexDecomposition_2`

*/
template <typename Kernel, typename Container>
class Small_side_angle_bisector_decomposition_2 {
public:

  /// @{
  /// @}

}; /* end Small_side_angle_bisector_decomposition_2 */
} /* end namespace CGAL */
