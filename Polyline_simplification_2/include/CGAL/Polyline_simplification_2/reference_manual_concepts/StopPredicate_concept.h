/// Models of this concept are passed to the polyline simplification algorithm to indicate
/// when to stop the process.
///
/// @heading Has Models:
/// - CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold
/// - CGAL::Polyline_simplification_2::Stop_below_count_threshold
/// - CGAL::Polyline_simplification_2::Stop_below_cost_threshold
class PolylineSimplificationStopPredicate
{
public :

  /// Indicates if the simplification must be stoped.
  /// This is called right before each vertex is about to be removed
  /// @param cdt The underlying constrained Delaunay triangulation which embeeds the polyline set
  /// @param q The current vertex about to be removed (given as a PolylineSimplificationVertex)
  /// @param cost The associated cost for removing the current vertex (as given by 'PolylineSimplificationCostFunction')
  /// @param initial_count The initial number of vertices in the entire polyline set (including intersection vertices not in any source polyline)
  /// @param current_count The current number of vertices
  /// @return true if the algorithm should stop, false if it should continue
  template<class ConstrainedDelaunayTriangulation, class VertexHandle>
  bool operator()( ConstrainedDelaunayTriangulation const& cdt
                 , VertexHandle const& q
                 , double cost
                 , std::size_t initial_count
                 , std::size_t current_count
                 ) const ;
};
