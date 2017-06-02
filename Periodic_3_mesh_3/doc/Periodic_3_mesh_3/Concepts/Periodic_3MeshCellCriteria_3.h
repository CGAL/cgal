/*!
\ingroup PkgPeriodic_3_mesh_3Concepts
\cgalConcept

\cgalRefines `MeshCellCriteria_3`

Similarly to the (non-periodic) 3D mesh generator, the Delaunay refinement
process involved in the template functions `make_periodic_3_mesh_3()`
and `refine_periodic_3_mesh_3()` is guided by a set of elementary refinement
criteria that concern either mesh tetrahedra or surface facets.
The concept `Periodic_3MeshCellCriteria_3` describes the types that
handle the refinement criteria for mesh tetrahedra.
From a syntaxic point of view, the concept `Periodic_3MeshCellCriteria_3` imposes
the same requirements as the concept `MeshCellCriteria_3`.

However, since periodic meshes are constructed within a single
fundamental domain, the oracles provided by this concept must be more powerful
than in `MeshCellCriteria_3` and handle periodicity.
For instance, the vertices of a cell that intersects the boundary
of the fundamental domain all live within the fundamental domain.
Directly using these coordinates to evaluate the badness of a facet would yield
incorrect results.
The badness oracle must thus be able to use translated images of vertices to
correctly test if a given cell satisfies a given criteria.

\cgalHasModel `CGAL::Periodic_3_mesh_cell_criteria_3<Tr>`

\sa `Periodic_3MeshCriteria_3`
\sa `CGAL::make_periodic_3_mesh_3()`
\sa `CGAL::refine_periodic_3_mesh_3()`

\sa `MeshCellCriteria_3`
*/

class Periodic_3MeshCellCriteria_3 {

/// @}

}; /* end Periodic_3MeshCellCriteria_3 */
