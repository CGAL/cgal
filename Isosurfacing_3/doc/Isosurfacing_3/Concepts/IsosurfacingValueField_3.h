/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

The concept `IsosurfacingValueField_3` describes the set of requirements to be fulfilled
by the value field template parameter of the domain classes `CGAL::Isosurfacing::Marching_cubes_domain_3` and `CGAL::Isosurfacing::Dual_contouring_domain_3`.

Value fields must be continuous and defined over the geometric span of the
space partitioning data structure (also known as "partition") being used.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Isosurfacing::Value_function_3}
\cgalHasModels{CGAL::Isosurfacing::Interpolated_discrete_values_3}
\cgalHasModelsEnd

\sa `IsosurfacingTraits_3`
\sa `IsosurfacingGradientField_3`
*/
class IsosurfacingValueField_3
{
public:
  /*!
  * The scalar type.
  */
  typedef unspecified_type FT;

  /*!
  * The 3D point type.
  */
  typedef unspecified_type Point_3;

  /*!
  * A descriptor that uniquely identifies a vertex (see `IsosurfacingPartition_3`).
  */
  typedef unspecified_type Vertex_descriptor;

  /*!
  returns the value of the value field at the point `p`.
  */
  FT operator()(const Point_3& p);

  /*!
  returns the value of the value field at the vertex `v`.
  */
  FT operator()(const Vertex_descriptor& v);
};
