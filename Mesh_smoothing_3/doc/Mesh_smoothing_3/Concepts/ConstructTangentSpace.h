/*!
\ingroup pkgMeshSmoothing3Concepts
\cgalConcept

The concept `ConstructTangentSpace` describes projection of mesh facets on target surfaces and 
mesh edge on target feature curves.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Mesh_smoothing_3::C3t3_mesh_projector}
\cgalHasModels{CGAL::Mesh_smoothing_3::C3t3_no_projection}
\cgalHasModelsEnd

\sa `CGAL::boundary_aware_mesh_smoothing`

*/
class ConstructTangentSpace {
private:

/// \name Types
/// @{

/*!
Model of `MeshComplex_3InTriangulation_3`
*/
using C3t3 = unspecified_type;

/*!
Geom_traits
*/
using Geom_traits = typename C3t3::Triangulation::Geom_traits;


/*!
Point type.
*/
using Point_3 = Geom_traits::Point_3;


/*!
Face associated with a patch
*/
using Patch_face = std::pair<C3t3::Surface_patch_index, C3t3::Facet>;

/*!
Edge associated with a curve
*/
using Curve_edge = std::pair<C3t3::Curve_index, C3t3::Edge>;


/// @}


/// \name Operations
/// The following functions are used to project entities
/// @{

/*!
returns a plane tangent to the patch to which its facet should align too. 
The list of Point_3 contains its current vertices location.  
*/
Tangent_space<Geom_traits> patch_face_projection_plane(Patch_face patch_face, std::vector<Point_3> face_points) const;


/*!
returns a line tangent to the curve to which its edge should align too. 
The array contains the two locations of its vertices.  
*/
Tangent_space<Geom_traits> curve_edge_projection_line(Curve_edge curve_edge, std::array<Point_3,2> edge_points) const;


/// @}



}; /* end ConstructTangentSpace */
