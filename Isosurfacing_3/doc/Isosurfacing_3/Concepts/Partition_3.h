/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

The concept `Partition_3` describes the set of requirements to be fulfilled
by the partition template parameter of the domain classes `CGAL::Isosurfacing::Marching_cubes_domain_3`
and `CGAL::Isosurfacing::Dual_contouring_domain_3`.

A 3D partition is a space partitioning data structure that provides a discrete representation
of a subset of 3D space.

A partial specialization of `CGAL::Isosurfacing::partition_traits` must be provided for all models.
This is similar to graph traits in \ref PkgBGL.

\cgalHasModelsBegin
\cgalHasModels{`CGAL::Isosurfacing::Cartesian_grid_3`}
\cgalHasModelsEnd
*/
class Partition_3
{
public:
  /*!
   * The geometric traits type.
   * Must be a model of `IsosurfacingTraits_3`.
  */
  typedef unspecified_type Geom_traits;

  /*!
   * \returns the geometric traits.
  */
  Geom_traits geom_traits();
};