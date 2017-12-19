/*!
\ingroup PkgPeriodic_3_mesh_3Concepts
\cgalConcept

\cgalRefines `MeshDomain_3`

The concept `Periodic_3MeshDomain_3` describes the knowledge required on the
object to be discretized.
The concept `Periodic_3MeshDomain_3` is the concept to be used when the input
domain is defined over the three-dimensional flat torus.

From a syntactic point of view, it defines almost exactly the same requirements
as the concept `MeshDomain_3` and thus `Periodic_3MeshDomain_3` refines `MeshDomain_3`.
However, the oracle must take into account the periodicity of the domain and
evaluate queries at the canonical representative a point (see Section
\ref Periodic_3_mesh_3InputDomain).

\cgalHasModel `CGAL::Implicit_periodic_3_mesh_domain_3<Function,BGT>`
\cgalHasModel `CGAL::Labeled_periodic_3_mesh_domain_3<LabelingFunction,BGT>`

\sa `CGAL::make_periodic_3_mesh_3()`
\sa `CGAL::refine_periodic_3_mesh_3()`
*/

class Periodic_3MeshDomain_3
{
public:
  /*!
  Cuboid type.
  */
  typedef unspecified_type Iso_cuboid_3;

  /*!
  Return the user-chosen cube that is the canonical instance of the flat torus.
  */
  const Iso_cuboid_3& canonical_periodic_domain();

/// @}

}; /* end Periodic_3MeshDomain_3 */
