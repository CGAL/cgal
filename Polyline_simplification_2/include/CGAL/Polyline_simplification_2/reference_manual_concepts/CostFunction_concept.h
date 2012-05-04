/// Models of this concept are passed to the polyline simplification algorithm to calculate
/// the "cost" of removing a vertex. Such a cost represents some measure of the deviation error between the
/// polyline sets before and after removal. The smaller the error the lower the cost. The algoritm processes
/// vertices in increasing cost order to preserve the overall polyline set shape as much as possible
///
/// @heading Has Models:
/// - CGAL::Polyline_simplification_2::Hybrid_squared_distance_cost
/// - CGAL::Polyline_simplification_2::Scaled_squared_distance_cost
/// - CGAL::Polyline_simplification_2::Squared_distance_cost
class PolylineSimplificationCostFunction
{

public:

  /// Given three consecutive polyline vertices "p,q,r", calculates the cost of removing
  /// vertex "q", replacing edges "p-q" and "q-r" with edge "p-r".
  /// Vertices "p,q,r" are consecutive in the current polyline, "possibly already simplified".
  /// @param cdt The underlying constrained Delaunay triangulation which embeeds the polyline set
  /// @param p The previous polyline vertex (given as a 'PolylineSimplificationVertex')
  /// @param q The current polyline vertex about to be removed
  /// @param r The next polyline vertex
  /// @param original_subpolyline_points_begin  Start iterator for the range of points that describe the original (unmodified) polyline subsequence from "p" to "r" (inclusive)
  /// @param original_subpolyline_points_end    Past-the-end iterator for the range of points that describe the original (unmodified) polyline subsequence from "p" to "r" (inclusive)
  /// @return The cost for removing "q". A result of "boost::none" can be used to indicate an infinite or uncomputable cost.
  template<class ConstrainedDelaunayTriangulation, class VertexHandle, class PointIterator>
  boost::optional<double> operator()( ConstrainedDelaunayTriangulation const& cdt
                                    , VertexHandle const& p
                                    , VertexHandle const& q
                                    , VertexHandle const& r
                                    , PointIterator           original_subpolyline_points_begin
                                    , PointIterator           original_subpolyline_points_end
                                    ) const ;

};

