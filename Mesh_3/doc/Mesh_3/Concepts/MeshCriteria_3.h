/*!
\ingroup PkgMesh3Concepts
\cgalConcept

The Delaunay refinement process involved in the
template functions `CGAL::make_mesh_3()` and `CGAL::refine_mesh_3()`
is guided by a set of elementary refinement criteria
that concern either mesh tetrahedra or surface facets.
The refinement criteria for tetrahedra are described
through the concept `MeshCellCriteria_3`
while the refinement criteria for surface facets
are described by the concept `MeshFacetCriteria_3`.
The concept `MeshCriteria_3` encapsulates these concepts.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Mesh_criteria_3<Tr>}
\cgalHasModelsEnd

\sa `MeshEdgeCriteria_3`
\sa `MeshFacetCriteria_3`
\sa `MeshCellCriteria_3`
\sa `MeshCriteriaWithFeatures_3`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`


*/

class MeshCriteria_3 {
public:

/// \name Types
/// @{

/*!
Functor that describes the criteria for
surface facets. This type must be a model
of the concept `MeshFacetCriteria_3`.
*/
typedef unspecified_type Facet_criteria;

/*!
Functor that describes the criteria for
mesh tetrahedra. This type must be a model of the concept
`MeshCellCriteria_3`.
*/
typedef unspecified_type Cell_criteria;

/// @}

/// \name Access Functions
/// @{

/*!
Returns the facet criteria.
*/
Facet_criteria facet_criteria_object();

/*!
Returns the cell criteria.
*/
Cell_criteria cell_criteria_object();

/// @}

}; /* end MeshCriteria_3 */
