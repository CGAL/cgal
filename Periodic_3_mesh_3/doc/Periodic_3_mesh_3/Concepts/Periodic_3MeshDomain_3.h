/*!
\ingroup PkgPeriodic_3_mesh_3Concepts
\cgalConcept

\cgalRefines `MeshDomain_3`

The concept `Periodic_3MeshDomain_3` describes the knowledge required on the
object to be discretized.
The concept `Periodic_3MeshDomain_3` is the concept to be used when the input
domain is defined over the three-dimensional flat torus. From a syntaxic point
of view, it defines exactly the same requirement as the concept `MeshDomain_3`.

However, since periodic meshes are constructed by considering a single fundamental
domain, the oracles must be more powerful than in `MeshDomain_3`
and handle periodicity.
For instance, when evaluating the `Do_intersect_surface` oracle for a segment
that intersects a predicate facet of the fundamental domain, it may be the case
that does not intersect the surface in the domain, but its translated copy that
intersects the opposite domain facet does intersect the surface in the domain.
The oracle must thus use translated images of such primitives to detect
if an intersection exists.

\cgalHasModel `CGAL::Implicit_periodic_3_mesh_domain_3<Function,BGT>`
\cgalHasModel `CGAL::Labeled_periodic_3_mesh_domain_3<LabelingFunction,BGT>`

\sa `CGAL::make_periodic_3_mesh_3()`
\sa `CGAL::refine_periodic_3_mesh_3()`
*/

class Periodic_3MeshDomain_3 {

/// @}

}; /* end Periodic_3MeshDomain_3 */
