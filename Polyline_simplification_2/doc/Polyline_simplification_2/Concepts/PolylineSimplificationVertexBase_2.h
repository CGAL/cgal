
/*!
\ingroup PkgPolylineSimplification2Concepts
\cgalConcept

The polyline simplification algorithm stores in the vertices
whether a vertex can be removed, and the cost of the removal.

\cgalRefines `TriangulationVertexBase_2` 

### Types ###

Defines the same types as the `TriangulationVertexBase_2` concept 

\cgalHasModel `CGAL::Polyline_simplification_2::Vertex_base_2<Vb>` 

\sa `TriangulationFaceBase_2` 
\sa `CGAL::Constrained_triangulation_plus_2<Tr>`

*/

class PolylineSimplificationVertexBase_2 {
public:
  /// A number type which must be the same as the <code>FT</code> of the geometric traits class of the triangulation
  typedef unspecified_type FT;

/// \name Access Functions 
/// @{

/*! indicates whether the vertex can be removed.
*/ 
bool is_removable() const; 

/*! allows to set whether the vertex can be removed.
*/
void set_removable(bool b);

/*! returns the cost of the vertex removal.
*/ 
FT cost() const; 

/*! allows to set the cost of the vertex removal.
*/
  void set_cost(const FT& ft);
/// @} 


}; /* end PolylineSimplificationVertexBase_2 */

