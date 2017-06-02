/*!
\ingroup PkgPeriodic_3_mesh_3Concepts
\cgalConcept

Similarly to the (non-periodic) 3D mesh generator, the Delaunay refinement
process involved in the template functions `make_periodic_3_mesh_3()`
and `refine_periodic_3_mesh_3()` is guided by a set of elementary refinement
criteria that concern either mesh tetrahedra or surface facets.
The refinement criteria for tetrahedra are described
through the concept `Periodic_3MeshCellCriteria_3`
while the refinement criteria for surface facets
are described by the concept `Periodic_3MeshFacetCriteria_3`.
The concept `Periodic_3MeshCriteria_3` encapsulates these concepts.

\cgalHasModel `CGAL::Periodic_3_mesh_criteria_3<Tr>`

\sa `MeshFacetCriteria_3`
\sa `MeshCellCriteria_3`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`
\sa `MeshCriteriaWithFeatures_3`

*/
class Periodic_3MeshCriteria_3 {

/// @}

}; /* end MeshCriteria_3 */
