
namespace CGAL {
namespace Surface_mesh_simplification {
/*!
\ingroup PkgSurfaceMeshSimplificationRef

The class `Surface_mesh_simplification::Edge_collapse_visitor_base` provides a base class for models
of the `EdgeCollapseSimplificationVisitor` concept.

This base class implements all of the visitor's callbacks.
This way, users need only override the callbacks they are interested in.
The callbacks <I>are not virtual</I> because this is not a dynamically polymorphic base class
and the derived visitor will never be used polymorphically at runtime (is perfectly fine to override
and hide a non-virtual method in the context of the static polymorphism used in the simplification algorithm).

\tparam TriangleMesh is the type of surface mesh being simplified, and must be a model
of the `MutableFaceGraph` and `HalfedgeListGraph` concepts.

\cgalModels `EdgeCollapseSimplificationVisitor`

*/
template <typename TriangleMesh>
struct Edge_collapse_visitor_base {

}; /* end Surface_mesh_simplification::Edge_collapse_visitor_base */

  } // namespace Surface_Mesh_Simplification
} // namespace CGAL
