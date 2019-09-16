namespace CGAL
{

/**
 * \ingroup PkgTriangulation3TraitsClasses
 *
 * @class Robust_weighted_circumcenter_filtered_traits_3
 *
 * Upgrades the functors models of `Kernel::ConstructWeightedCircumcenter_3`,
 * `Kernel::ComputeSquaredRadius_3`, and `Kernel::ComputeSquaredRadiusSmallestOrthogonalSphere_3`
 * that are provided by `K` with robust versions. This improved robustness
 * is achieved by using exact computations if the geometric configuration
 * is close to degenerate (e.g. almost coplanar points). The upgrade is completely
 * invisible from an interface point of view as the class `Robust_weighted_circumcenter_filtered_traits_3`
 * overrides the types and function objects associated with the basic versions.
 *
 * \tparam K must be a model of the `Kernel` concept.
 *
 * \sa \ref devman_robustness
 */
template<class K>
class Robust_weighted_circumcenter_filtered_traits_3
  : public K
{
public:
  /// \name Types
  /// @{

  /*!
    a model of `Kernel::ConstructWeightedCircumcenter_3`
  */
  typedef unspecified_type Construct_weighted_circumcenter_3;

  /*!
    a model of `Kernel::ComputeSquaredRadius_3`
  */
  typedef unspecified_type Compute_squared_radius_3;

  /*!
    a model of `Kernel::ComputeSquaredRadiusSmallestOrthogonalSphere_3`
  */
  typedef unspecified_type Compute_squared_radius_smallest_orthogonal_sphere_3;

  /// @}

  /// \name Operations
  /// @{

  Construct_weighted_circumcenter_3
  construct_weighted_circumcenter_3_object();

  Compute_squared_radius_3
  compute_squared_radius_3_object();

  Compute_squared_radius_smallest_orthogonal_sphere_3
  compute_squared_radius_smallest_orthogonal_sphere_3_object();

  /// @}
};

}  // end namespace CGAL
