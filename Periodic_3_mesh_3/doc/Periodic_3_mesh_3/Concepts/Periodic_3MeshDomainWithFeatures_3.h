/*!
\ingroup PkgPeriodic_3_mesh_3Concepts
\cgalConcept

\cgalRefines `MeshDomainWithFeatures_3` and `Periodic_3MeshDomain_3`

The concept `Periodic_3MeshDomainWithFeatures_3` describes the knowledge required on the
object to be discretized.

While the concept `Periodic_3MeshDomain_3` only exposes the 2-dimensional and
3-dimensional features of the periodic domain through different queries,
the concept `Periodic_3MeshDomainWithFeatures_3` also exposes 0 and 1-dimensional features.
The exposed features of the domain are respectively called subdomains, surface patches,
 curves and corners according to their respective dimensions 3, 2, 1, and 0.

From a syntactic point of view, it defines exactly the same requirement
as the concept `MeshDomainWithFeatures_3` and thus `Periodic_3MeshDomainWithFeatures_3`
refines `MeshDomainWithFeatures_3` without any additional requirement.
However, the oracle must take into account the periodicity of the domain.

The class `CGAL::Mesh_domain_with_polyline_features_3` is a model of the concept
`MeshDomainWithFeatures_3` and thus wrapping any model of `Periodic_3MeshDomain`
with `CGAL::Mesh_domain_with_polyline_features_3` gives a model of `Periodic_3MeshDomainWithFeatures_3`.

\cgalHasModel `CGAL::Mesh_domain_with_polyline_features_3<
                 CGAL::Implicit_periodic_3_mesh_domain_3<Function,BGT> >`
\cgalHasModel `CGAL::Mesh_domain_with_polyline_features_3<
                 CGAL::Labeled_periodic_3_mesh_domain_3<LabelingFunction,BGT> >`

\sa `CGAL::make_periodic_3_mesh_3()`
\sa `CGAL::refine_periodic_3_mesh_3()`
*/

class Periodic_3MeshDomainWithFeatures_3 {

/// @}

}; /* end Periodic_3MeshDomainWithFeatures_3 */
