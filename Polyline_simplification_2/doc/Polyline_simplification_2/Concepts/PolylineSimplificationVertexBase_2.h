
/*!
\ingroup PkgPolylineSimplification2Concepts
\cgalConcept

Thr polyline simplification algorithm stores in the vertices
whether a vertex can be removed, and the cost of the removal.

\cgalRefines `TriangulationVertexBase_2` 

### Types ###

Defines the same types as the `TriangulationVertexBase_2` concept 

\cgalHasModel `CGAL::Polyline_simplification_2::Vertex_base_2<Vb>` 

\sa `TriangulationFaceBase_2` 
\sa `CGAL::Polyline_constrained_triangulation_2<CDT>` 

*/

class PolylineSimplificationVertexBase_2 {
public:
  /// A number type.
  typedef Hidden_type FT;

/// \name Access Functions 
/// @{

/*! The Boolean that indicates whether the vertex can be removed.
*/ 
bool& unremovable(); 

/*! The cost if the vertex got removed.
*/ 
FT& cost(); 

/// @} 


}; /* end PolylineSimplificationVertexBase_2 */

