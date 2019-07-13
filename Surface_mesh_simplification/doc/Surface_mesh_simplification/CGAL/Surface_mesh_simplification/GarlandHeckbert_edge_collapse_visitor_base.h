namespace CGAL {
namespace Surface_mesh_simplification {
/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `Surface_mesh_simplification::GarlandHeckbert_edge_collapse_visitor_base`
provides a base class for models of the `EdgeCollapseSimplificationVisitor` concept
to be used as a visitor for Garland-Heckbert strategy.

Derived classes should call the contructor of this visitor with an externally provided
state, and a discontinuity multiplier. Derived classes should also call OnStarted, and
OnCollapsed functions of this class if they choose to override them.

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.

\cgalModels `EdgeCollapseSimplificationVisitor`

\sa `Surface_mesh_simplification::GarlandHeckbert_cost<TriangleMesh>`

\sa `Surface_mesh_simplification::GarlandHeckbert_placement<TriangleMesh>`

*/
template<class TriangleMesh>
struct GarlandHeckbert_edge_collapse_visitor_base : Edge_collapse_visitor_base<TriangleMesh> {
  typedef unspecified_type garland_heckbert_state_type;
  typedef unspecified_type FT;

  /// \name Creation
  /// @{

  /*!
  Initilizes the visitor with given the <I>garland heckbert state</I> object, and
  the <I>discontinuity multiplier</I>.
  Garland&Heckbert strategy requires a shared state object between cost, placement, and visitor policies.
  Discontinuity multiplier is used to increase costs at borders of meshes. This way, collapsing vertices on borders
  gets penalized heavily.

  */
  GarlandHeckbert_edge_collapse_visitor_base(
    garland_heckbert_state_type& state, FT discontinuity_multiplier = 100.0);
  /// @}
};
}
}
