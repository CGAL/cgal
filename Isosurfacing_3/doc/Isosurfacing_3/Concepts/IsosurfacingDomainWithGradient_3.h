/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

\cgalRefines{IsosurfacingDomain_3}

\brief The concept `IsosurfacingDomainWithGradient_3` describes the set of requirements to be
fulfilled by any class used as input data for some isosurfacing algorithms.

This concept refines `IsosurfacingDomain_3` to add a `gradient()` function which is used
by isosurfacing domains to query the domain for the gradient of the values field
at a 3D query point (not necessarily a vertex) in space.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Isosurfacing::Dual_contouring_domain_3}
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
  The 3D point type.
  */
  typedef Geom_traits::Point_3 Point_3;

  /*!
  The 3D vector type.
  */
  typedef Geom_traits::Vector_3 Vector_3;

  /// @}

  /// \name Operations
  /// @{

  /*!
  returns the gradient at the position `p`
  */
  Vector_3 gradient(const Point_3& p) const;

  /// @}
};
