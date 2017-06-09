/*!
\ingroup PkgPeriodic_3_mesh_3Concepts
\cgalConcept

The function object concept `Periodic_3MeshEdgeCriteria_3` is designed to drive
the process which samples the 1-dimensional features of the domain.
It provides an upper bound for the distance between two protecting ball centers
that are consecutive on a 1-feature.

From a syntaxic point of view, the concept `Periodic_3MeshEdgeCriteria_3` imposes
the same requirements as the concept `MeshEdgeCriteria_3`.
However, since periodic meshes are constructed within a single
fundamental domain, the oracles provided by this concept must be more powerful
than in `MeshEdgeCriteria_3` and handle periodicity.
For instance, the vertices of an edge that intersects the boundary
of the fundamental domain all live within the fundamental domain.
Directly using these coordinates to evaluate the length of an edge would yield
incorrect results.
The oracle must thus use be able to use translated images of vertices to
correctly evaluate the length of edges.
\cgalHasModel `CGAL::Periodic_3_mesh_edge_criteria_3<Tr>`

\sa `Periodic_3MeshCriteria_3`

\sa `MeshEdgeCriteria_3`
\sa `MeshCriteriaWithFeatures_3`
*/

class Periodic_3MeshEdgeCriteria_3 {

/// @}

}; /* end Periodic_3MeshEdgeCriteria_3 */
