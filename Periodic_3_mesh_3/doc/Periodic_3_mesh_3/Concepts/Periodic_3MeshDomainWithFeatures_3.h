/*!
\ingroup PkgPeriodic3Mesh3Concepts
\cgalConcept

\cgalRefines `MeshDomainWithFeatures_3` and `Periodic_3MeshDomain_3`

The concept `Periodic_3MeshDomainWithFeatures_3` describes the knowledge required on the
object to be discretized.

While the concept `Periodic_3MeshDomain_3` only exposes the 2-dimensional and
3-dimensional features of the periodic domain through different queries,
the concept `Periodic_3MeshDomainWithFeatures_3` also exposes 0 and 1-dimensional features.
The exposed features of the domain are respectively called subdomains, surface patches,
curves, and corners according to their respective dimensions 3, 2, 1, and 0.

From a syntactic point of view, `Periodic_3MeshDomainWithFeatures_3`
refines `MeshDomainWithFeatures_3`. However, the various requirements from
`MeshDomainWithFeatures_3` must also take into account the periodicity of the domain
(see Section \ref Periodic_3_mesh_3InputDomain).

Wrapping any model of `Periodic_3MeshDomain_3` with the class
`CGAL::Mesh_domain_with_polyline_features_3` gives a model
of `Periodic_3MeshDomainWithFeatures_3`.

\cgalHasModel `CGAL::Mesh_domain_with_polyline_features_3<
                 CGAL::Labeled_mesh_domain_3<BGT> >`

\sa `CGAL::Periodic_3_function_wrapper<Function,BGT>`

\sa `CGAL::make_periodic_3_mesh_3()`
\sa `CGAL::refine_periodic_3_mesh_3()`
*/

class Periodic_3MeshDomainWithFeatures_3
{
public:
  /*!
  Returns the indices of the curves incident to the corner `id`, if any.

  \tparam IndicesOutputIterator must meet the requirements of `OutputIterator`,
                                with value type `MeshDomainWithFeatures_3::Curve_index`.
  */
  template <typename IndicesOutputIterator>
  IndicesOutputIterator
  get_corner_incident_curves(Corner_index id, IndicesOutputIterator out);

}; /* end Periodic_3MeshDomainWithFeatures_3 */
