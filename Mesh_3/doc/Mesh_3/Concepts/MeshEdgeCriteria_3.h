/*!
\ingroup PkgMesh3Concepts
\cgalConcept

The function object concept `MeshEdgeCriteria_3` is designed to drive the process which samples
the 1-dimensional features of the domain.
It provides an upper bound for the distance between two protecting ball centers
that are consecutive on a 1-feature.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Mesh_edge_criteria_3<Tr>}
\cgalHasModelsEnd

\sa `MeshCellCriteria_3`
\sa `MeshFacetCriteria_3`
\sa `MeshCriteria_3`
\sa `MeshCriteriaWithFeatures_3`

*/

class MeshEdgeCriteria_3 {
public:

/// \name Types
/// @{

/*!
Point type. Must match the `Point_3` type in
the triangulation type used by the mesh generation function.
*/
typedef unspecified_type Point_3;

/*!
Type for edges of the triangulation. Must match the
`Edge` type in the triangulation type used by the mesh generation function.
*/
typedef unspecified_type Edge;

/*!
Numerical type.
*/
typedef unspecified_type FT;

/*!
Feature index type.
*/
typedef unspecified_type Index;

/// @}

/// \name Operations
/// @{

/*!

Returns the value of the sizing field (i.e., the maximum edge length) at point `p`,
lying on subcomplex of dimension `dim` and index `index`.
*/
FT sizing_field(const Point_3& p, const int dim, const Index& index);

/*!
Returns the lower bound on edge length, set by the user.
The lower bound is ignored when its value is 0.
*/
const FT& min_length_bound() const;

/*!
Returns the value of the distance field (i.e., the maximum edge distance) at point `p`
lying on subcomplex of dimension `dim` and index `index`.
*/
FT distance_field(const Point_3& p, const int dim, const Index& index);

/*!
Returns whether or not the distance field should be checked during the protection phase.
If false, the distance field is ignored.
*/
bool has_distance_field() const;
/// @}

}; /* end MeshEdgeCriteria_3 */
