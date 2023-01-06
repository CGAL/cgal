/*!
\ingroup PkgIsosurfacing3Concepts

\cgalConcept

The concept `IsosurfacingDomainWithGradient` describes the set of requirements to be
fulfilled by any class used as input data for some isosurfacing algorithms.

\cgalHasModel `CGAL::Isosurfacing::Explicit_cartesian_grid_domain`
\cgalHasModel `CGAL::Isosurfacing::Implicit_cartesian_grid_domain`
*/
class IsosurfacingDomainWithGradient
  : public IsosurfacingDomain
{
public:
  /// \name Types
  /// @{

  /*!
  The vector type.
  */
  typedef unspecified_type Vector;

  /// @}

  /// \name Operations
  /// The following member function must exist.
  /// @{

  /*!
  Get the gradient at a position

  \param p the point at which the gradient is evaluated

  \return the gradient vector
  */
  Vector gradient(const Point& p) const;

  /// @}
};
