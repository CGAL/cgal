/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

\cgalRefines `IsosurfacingDomain_3`

\brief The concept `IsosurfacingDomainWithGradient_3` describes the set of requirements to be
fulfilled by any class used as input data for some isosurfacing algorithms.

\cgalHasModel `CGAL::Isosurfacing::Explicit_Cartesian_grid_domain_3`
\cgalHasModel `CGAL::Isosurfacing::Implicit_Cartesian_grid_domain_3`
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
  The vector type.
  */
  typedef Geom_traits::Vector_3 Vector_3;

  /// @}

  /// \name Operations
  /// @{

  /*!
  gets the gradient at a position
  */
  Vector_3 gradient(const Geom_traits::Point_3& p) const;

  /// @}
};
