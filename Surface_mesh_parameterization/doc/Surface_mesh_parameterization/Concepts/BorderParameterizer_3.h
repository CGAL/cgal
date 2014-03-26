
/*!
\ingroup PkgSurfaceParameterizationConcepts
\cgalConcept

`BorderParameterizer_3` is a concept of class that parameterizes a given type of mesh, `Adaptor`, which is a model of the `ParameterizationMesh_3` concept. 

Implementation note: To simplify the implementation, `BorderParameterizer_3` models know only the `ParameterizationMesh_3` class. They do not know the parameterization algorithm requirements or the kind of sparse linear system used. 

Design Pattern 
-------------- 

`BorderParameterizer_3` models are Strategies \cgalCite{cgal:ghjv-dpero-95} : they implement a strategy of border parameterization for models of `ParameterizationMesh_3`. 

Creation 
-------------- 

Construction and destruction are undefined. 

\cgalHasModel `CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>`
\cgalHasModel `CGAL::Circular_border_uniform_parameterizer_3<ParameterizationMesh_3>`
\cgalHasModel `CGAL::Square_border_arc_length_parameterizer_3<ParameterizationMesh_3>` 
\cgalHasModel `CGAL::Square_border_uniform_parameterizer_3<ParameterizationMesh_3>`
\cgalHasModel `CGAL::Two_vertices_parameterizer_3<ParameterizationMesh_3>`

\sa `ParameterizerTraits_3`
\sa `ParameterizationMesh_3`

*/

class BorderParameterizer_3 {
public:

/// \name Types 
/// @{

/*!

Export `ParameterizationMesh_3` template parameter. 

*/ 
typedef unspecified_type Adaptor; 

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
Error_code parameterize_border(Adaptor& mesh); 

/*!

Indicate if border's shape is convex. 

*/ 
bool is_border_convex(); 

/// @}

}; /* end BorderParameterizer_3 */

