
/*!
\ingroup PkgSurfaceParameterizationConcepts
\cgalConcept

`BorderParameterizer_3` is a concept of class that parameterizes the border of
a given type of mesh, `TriangleMesh`, which is a model of the `FaceGraph` concept.
Some of the vertices (possibly all) of the mesh are assigned a 2D position on a
shape (a geometrical object).

\cgalHasModel `CGAL::Circular_border_parameterizer_3<TriangleMesh>`
\cgalHasModel `CGAL::Square_border_parameterizer_3<TriangleMesh>`
\cgalHasModel `CGAL::Two_vertices_parameterizer_3<TriangleMesh>`

\sa `ParameterizerTraits_3`

*/

class BorderParameterizer_3 {
public:

/// \name Types
/// @{
/*!

A given polygon mesh type, TriangleMesh, which is a model of the `FaceGraph` concept.

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

Assign a 2D position (i.e.\ a `(u, v)` pair) on the shape  to (some of) the vertices on the
border of the mesh. Mark them as <I>parameterized</I>.

*/
Error_code parameterize_border(const TriangleMesh& mesh);

/*!

Indicate if the shape is convex.

*/
bool is_border_convex() const;

/// @}

}; /* end BorderParameterizer_3 */

