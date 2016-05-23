
/*!
\ingroup PkgSurfaceParameterizationConcepts
\cgalConcept

`BorderParameterizer_3` is a concept of class that parameterizes a given type of mesh, `TriangleMesh`, which is a model of the `FaceGraph` concept. 


\cgalHasModel `CGAL::Circular_border_arc_length_parameterizer_3<TriangleMesh>`
\cgalHasModel `CGAL::Circular_border_uniform_parameterizer_3<TriangleMesh>`
\cgalHasModel `CGAL::Square_border_arc_length_parameterizer_3<TriangleMesh>` 
\cgalHasModel `CGAL::Square_border_uniform_parameterizer_3<TriangleMesh>`
\cgalHasModel `CGAL::Two_vertices_parameterizer_3<TriangleMesh>`

\sa `ParameterizerTraits_3`

*/

class BorderParameterizer_3 {
public:

/// \name Types 
/// @{
/*!


*/ 
  typedef unspecified_type TriangleMesh; 
/*!

The various errors detected by this package. 

*/ 
typedef unspecified_type Error_code; 

/// @} 

/// \name Operations 
/// @{

/*!

Assign to mesh's border vertices a 2D position (i.e.\ a `(u, v)` pair) on border's shape. Mark them as <I>parameterized</I>. Return false on error. 

*/ 
Error_code parameterize_border(const TriangleMesh& mesh); 

/*!

Indicate if border's shape is convex. 

*/ 
bool is_border_convex(); 

/// @}

}; /* end BorderParameterizer_3 */

