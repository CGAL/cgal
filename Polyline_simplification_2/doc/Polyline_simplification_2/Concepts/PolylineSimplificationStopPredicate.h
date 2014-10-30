/*!
\ingroup PkgPolylineSimplification2Concepts
\cgalConcept

Models of this concept are passed to the polyline simplification algorithm to indicate
when to stop the process.

\cgalHasModel `CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold`
\cgalHasModel `CGAL::Polyline_simplification_2::Stop_below_count_threshold`
\cgalHasModel `CGAL::Polyline_simplification_2::Stop_above_cost_threshold`
*/

class PolylineSimplificationStopPredicate
{
public :
  /*!
  Indicates if the simplification must be stopped.
  This is called right before each vertex is about to be removed
  \param ct The underlying constrained Delaunay triangulation with constraint hierarchy which embeds the polyline constraints
  \param q The current vertex about to be removed
  \param cost The associated cost for removing the current vertex (as given by `PolylineSimplificationCostFunction`)
  \param initial_count The initial number of vertices in the entire polyline set (including intersection vertices not in any source polyline)
  \param current_count The current number of vertices
  \return `true` if the algorithm should stop, `false` if it should continue.
\tparam Tr must be `CGAL::Constrained_Delaunay_triangulation_2` with a vertex type that
is model of `PolylineSimplificationVertexBase_2`.
  */
  template<class Tr>
  bool operator()( const CGAL::Constrained_triangulation_plus_2<Tr>& ct
                  , CGAL::Constrained_triangulation_plus_2<Tr>::Vertex_handle q
                 , double cost
                 , std::size_t initial_count
                 , std::size_t current_count
                 ) const ;
};
