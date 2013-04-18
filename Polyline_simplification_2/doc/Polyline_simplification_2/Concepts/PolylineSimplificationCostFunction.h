
/*!
\ingroup PkgPolylineSimplification2Concepts
\cgalConcept

Models of this concept are passed to the polyline simplification
algorithm to calculate the *cost* of removing a vertex. Such a cost
represents some measure of the deviation error between the polyline
sets before and after removal. The smaller the error the lower the
cost. The algorithm processes vertices in increasing cost order to
preserve the overall polyline set shape as much as possible


\cgalHasModel `CGAL::Polyline_simplification_2::Hybrid_squared_distance_cost`
\cgalHasModel `CGAL::Polyline_simplification_2::Scaled_squared_distance_cost`
\cgalHasModel `CGAL::Polyline_simplification_2::Squared_distance_cost`

*/

class PolylineSimplificationCostFunction {
public:


/*! 
Given three consecutive polyline vertices `*vip, *viq, *vir`, calculates the cost of removing vertex `*viq`, replacing edges `(*vip,*viq)` and `(*viq,*vir)` with edge `(*vip,*vir)`. 

\param pct The underlying polyline constrained  triangulation which embeds the polyline set.
\param vip The first vertex
\param viq The second vertex
\param vir The third vertex
\returns The cost for removing `*viq`. A result of `boost::none` can be used to indicate an infinite or uncomputable cost.

`Tr::Geom_traits` must provide a functor `Compute_squared_distance` with an operator `Tr::Geom_traits::FT operator()(Tr::Point, Tr::Point)`.
*/ 
boost::optional<Polyline_constrained_triangulation_2<Tr>::Geom_traits::FT> 
          operator()(Polyline_constrained_triangulation_2<Tr> const& pct, 
                     Polyline_constrained_triangulation_2<Tr>::Vertices_in_constraint_iterator vip, 
                     Polyline_constrained_triangulation_2<Tr>::Vertices_in_constraint_iterator viq, 
                     Polyline_constrained_triangulation_2<Tr>::Vertices_in_constraint_iterator vir) const;}

}; /* end TriangulationVertexBase_2 */

