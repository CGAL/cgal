
namespace CGAL {

/*!
\ingroup PkgSurfaceMeshSimplification

The class `Surface_mesh_simplification::Edge_collapse_visitor_base` provides a base class for models of the `EdgeCollapseSimplificationVisitor` concept. 
It has one template argument: the type of surface being simplified, 
which must be a model of the `EdgeCollapsableMesh` concept. 

This base class implements all of the visitor's callbacks. 
This way, users need only override the callbacks they are interested in. 
The callbacks <I>are not virtual</I> because this is not a dynamically polymorphic base class 
and the derived visitor will never be used polymorphically at runtime (is perfectly fine to override 
and hide a non-virtual method in the context of the static polymorphism used in the simplification algorithm). 

\models ::EdgeCollapseSimplificationVisitor 

*/
template< typename ECM >
class Surface_mesh_simplification::Edge_collapse_visitor_base {
public:

/// @}

}; /* end Surface_mesh_simplification::Edge_collapse_visitor_base */
} /* end namespace CGAL */
