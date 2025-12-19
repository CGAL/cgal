
/*!
\ingroup PkgPolylineSimplification2Concepts
\cgalConcept

Models of this concept are passed to the polyline simplification
algorithm to calculate the *cost* of removing a vertex. Such a cost
represents some measure of the deviation error between the polyline
sets before and after removal. The smaller the error the lower the
cost. The algorithm processes vertices in increasing cost order to
preserve the overall polyline set shape as much as possible

\cgalRefines{CopyConstructible,Assignable}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Polyline_simplification_2::Hybrid_squared_distance_cost}
\cgalHasModels{CGAL::Polyline_simplification_2::Scaled_squared_distance_cost}
\cgalHasModels{CGAL::Polyline_simplification_2::Squared_distance_cost}
\cgalHasModelsEnd

*/

class PolylineSimplificationCostFunction {
public:


/*!
Given a vertex in constraint iterator `viq` computes `vip=std::prev(viq)` and `vir=std::next(vir)`, and the cost of removing vertex `*viq`, replacing edges `(*vip,*viq)` and `(*viq,*vir)` with edge `(*vip,*vir)`.

\param ct The underlying constrained Delaunay triangulation which embeds the polyline constraints
\param viq The vertex in constraint iterator of the vertex to remove
\returns The cost for removing `*viq`. The value `std::nullopt` can be returned to indicate an infinite or uncomputable cost.

\tparam CDT must be `CGAL::Constrained_triangulation_plus_2` with a vertex type that
is model of `PolylineSimplificationVertexBase_2`. `CDT::Geom_traits` must be model of
the concept `ConstrainedDelaunayTriangulationTraits_2`.

*/
  template <typename CDT>
  std::optional<CDT::Geom_traits::FT>
  operator()(CGAL::Constrained_triangulation_plus_2<CDT> const& ct,
             CGAL::Constrained_triangulation_plus_2<CDT>::Vertices_in_constraint_iterator viq) const;}

}; /* end TriangulationVertexBase_2 */

