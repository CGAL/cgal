/*!
\ingroup pkgMeshSmoothing3Concepts
\cgalConcept

The concept `C3t3Projector` describes projection of mesh facets on target surfaces and 
mesh edge on target feature curves.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Mesh_smoothing_3::C3t3_mesh_projector}
\cgalHasModels{CGAL::Mesh_smoothing_3::C3t3_no_projection}
\cgalHasModelsEnd

\sa `CGAL::boundary_aware_mesh_smoothing`

*/
class C3t3Projector {
public:

/// \name Types
/// @{

/*!
Point type.
*/
using Point_3 = unspecified_type;


/*!
Vector type.
*/
using Vector_3 = unspecified_type;

/*!
Surface patch index
*/
using Surface_patch_index = unspecified_type;

/*!
Face descriptor
*/
using Facet = unspecified_type;

/*!
Face associated with a patch
*/
using Patch_face = std::pair<Surface_patch_index, Facet>;


/*!
Curve index
*/
using Curve_index = unspecified_type;

/*!
Edge descriptor
*/
using Edge = unspecified_type;

/*!
Edge associated with a curve
*/
using Curve_edge = `std::pair<Curve_index, Edge>`;


/// @}

/*!
Return the plane the patch face should align to. 
*/
`std::pair<Point_3, Vector_3>` patch_projection_plane(Patch_face patch_face, Point_3 face_center, double face_radius) const;

/*!
Return if a patch face should be projected or not. 
*/
bool project_patch_face(Patch_face patch_face) const;


/*!
Return the line the curve edge should align to. 
*/
std::pair<Point_3, Vector_3> curve_projection_tangent(Curve_edge curve_edge, Point_3 edge_center, double segment_size) const;

/*!
Return if a curve edge should be projected or not. 
*/
bool project_curve_edge(Curve_edge curve_edge) const,


/// @}



}; /* end MeshDataStructure */
