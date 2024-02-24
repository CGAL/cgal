/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

The concept `GradientField_3` describes the set of requirements to be fulfilled
by the gradient field template parameter of the domain class `CGAL::Isosurfacing::Dual_contouring_domain_3`.

Gradient fields must be continuous and defined over the geometric span of the
space partitioning data structure (also known as "partition") being used.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Isosurfacing::Gradient_function_3}
\cgalHasModels{CGAL::Isosurfacing::Finite_difference_gradient_3}
\cgalHasModels{CGAL::Isosurfacing::Interpolated_discrete_gradients_3}
\cgalHasModelsEnd

\sa `IsosurfacingTraits_3`
\sa `ValueField_3`
*/
class GradientField_3
{
public:
  /*!
  * The 3D point type.
  */
  typedef unspecified_type Point_3;

  /*!
  * The 3D vector type.
  */
  typedef unspecified_type Vector_3;

  /*!
  returns the gradient at the position `p`
  */
  Vector_3 operator()(const Point_3& p) const;
};
