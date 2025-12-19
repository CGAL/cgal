/*!
\ingroup PkgMesh3Concepts
\cgalConcept

The concept `MeshCriteriaWithFeatures_3` refines the concept `MeshCriteria_3`.
The concept `MeshCriteria_3` encapsulates
the concepts `MeshCellCriteria_3` and
`MeshFacetCriteria_3` describing the refinement
criteria for, respectively, mesh cells and surface facets.
For domains with features, the concept `MeshCriteriaWithFeatures_3`
additionally encapsulates the
concept `MeshEdgeCriteria_3`,
that describes the requirements, in terms of sizing, for the discretization of the domain \f$ 1\f$-dimensional features.

\cgalRefines{MeshCriteria_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Mesh_criteria_3<Tr>}
\cgalHasModelsEnd

\sa `MeshEdgeCriteria_3`
\sa `MeshFacetCriteria_3`
\sa `MeshCellCriteria_3`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/

class MeshCriteriaWithFeatures_3 {
public:

/// \name Types
/// @{

/*!
Functor that describes the criteria for
the mesh edges that discretize the input domain 1-dimensional features.
This type must be a model
of the concept `MeshEdgeCriteria_3`.
*/
typedef unspecified_type Edge_criteria;

/// @}

/// \name Access Functions
/// @{

/*!
Returns the edge criteria.
*/
Edge_criteria edge_criteria_object();

/// @}

}; /* end MeshCriteriaWithFeatures_3 */
