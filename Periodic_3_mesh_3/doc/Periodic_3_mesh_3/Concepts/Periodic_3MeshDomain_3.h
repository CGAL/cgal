/*!
\ingroup PkgPeriodic3Mesh3Concepts
\cgalConcept

\cgalRefines{MeshDomain_3}

The concept `Periodic_3MeshDomain_3` describes the knowledge required on the
object to be discretized.
The concept `Periodic_3MeshDomain_3` is the concept to be used when the input
domain is defined over the three-dimensional flat torus.

From a syntactic point of view, it defines almost the same requirements
as the concept `MeshDomain_3` and thus `Periodic_3MeshDomain_3` refines `MeshDomain_3`:
the concept `Periodic_3MeshDomain_3` additionally requires an access to the user-defined
canonical cuboid via the function `bounding_box`.
However, the oracle must take into account the periodicity of the domain (see Section
\ref Periodic_3_mesh_3InputDomain).

The class `CGAL::Labeled_mesh_domain_3<BGT>` paired with a periodic labeling function
is a model of this concept. It is possible to create artificially periodic functions
through the class `CGAL::Periodic_3_function_wrapper<Function,BGT>`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Labeled_mesh_domain_3<BGT>}
\cgalHasModelsEnd

\sa `CGAL::Labeled_mesh_domain_3<BG>`

\sa `CGAL::make_periodic_3_mesh_3()`
\sa `CGAL::refine_periodic_3_mesh_3()`
*/

class Periodic_3MeshDomain_3
{
public:
  /*!
  The canonical cuboid type.
  */
  typedef unspecified_type Iso_cuboid_3;

  /*!
  returns the user-chosen cuboid that is the canonical instance of the flat torus.
  */
  const Iso_cuboid_3& bounding_box();

}; /* end Periodic_3MeshDomain_3 */
