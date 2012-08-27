/*!
\ingroup PkgMesh_3Concepts
\cgalconcept

The function object concept `MeshEdgeCriteria_3` is designed to drive the process which samples 
the 1-dimensional features of the domain. 
It provides an upper bound for the distance between two protecting ball centers 
that are consecutive on a 1-feature. 

\hasModel `Mesh_edge_criteria_3<Tr>` 

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
typedef Hidden_type Point_3; 

/*! 
Type for edges of the triangulation. Must match the 
`Edge` type in the triangulation type used by the mesh generation function. 
*/ 
typedef Hidden_type Edge; 

/*! 
Numerical type. 
*/ 
typedef Hidden_type FT; 

/// @} 

/// \name Operations 
/// @{

/*! 

Returns `true` if edge `e` does not fulfill the criteria. 
*/ 
bool operator()(const Edge& e); 

/*! 

Returns the value of the sizing field (i.e. the maximum edge length) at point `p`. 
*/ 
FT sizing_field(const Point_3& p); 

/// @}

}; /* end MeshEdgeCriteria_3 */
