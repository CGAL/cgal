/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

\brief The concept `IsosurfacingDomainWithGradient_3` describes the set of requirements to be
fulfilled by any class used as input data for some isosurfacing algorithms.

This concept refines `IsosurfacingDomain_3` to add a `gradient()` function which is used
by isosurfacing domains to query the domain for the gradient of the implicit field
at a 3D query point (not necessarily a vertex) in space.

\cgalRefines `IsosurfacingDomain_3`

\cgalHasModelsBegin
\cgalHasModels{CGAL::Isosurfacing::Explicit_Cartesian_grid_domain_3}
\cgalHasModels{CGAL::Isosurfacing::Implicit_Cartesian_grid_domain_3}
\cgalHasModelsEnd
*/
class IsosurfacingDomainWithGradient_3
{
public:
  /// \name Types
  /// @{

  /*!
  The geometric traits type.
  Must be a model of `IsosurfacingTraits_3`.
  */
  typedef unspecified_type Geom_traits;

  /*!
  The point type.
  */
  typedef Geom_traits::Point_3 Point_3;

  /*!
  The vector type.
  */
  typedef Geom_traits::Vector_3 Vector_3;

  /// @}

  /// \name Operations
  /// @{

  /*!
  gets the gradient at the position `p`
  */
  Vector_3 gradient(const Point_3& p) const;

  /// @}
};
