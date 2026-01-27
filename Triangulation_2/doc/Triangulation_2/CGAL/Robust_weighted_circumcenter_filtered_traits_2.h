namespace CGAL
{

/**
 * \ingroup PkgTriangulation2TraitsClasses
 *
 * @class Robust_circumcenter_filtered_traits_2
 *
 * Upgrades the functors models of `Kernel::ConstructCircumcenter_2` and
 * `Kernel::ComputeSquaredRadius_2` that are provided by `K` with robust versions.
 * This improved robustness is achieved by using exact computations if the geometric configuration
 * is close to degenerate (e.g. almost coplanar points). The upgrade is completely
 * invisible from an interface point of view as the class `Robust_circumcenter_filtered_traits_2`
 * overrides the types and function objects provided by `K`.
 *
 * \tparam K must be a model of the `Kernel` concept.
 *
 * \sa `CGAL::Robust_weighted_circumcenter_filtered_traits_2<K>`
 * \sa \ref devman_robustness
 */
template<class K>
class Robust_circumcenter_filtered_traits_2
  : public K
{
public:
  /// \name Types
  /// @{

  /*!
    a model of `Kernel::ConstructCircumcenter_2`
  */
  typedef unspecified_type Construct_circumcenter_2;

  /*!
    a model of `Kernel::ComputeSquaredRadius_2`
  */
  typedef unspecified_type Compute_squared_radius_2;

  /// @}

  /// \name Operations
  /// @{

  Construct_circumcenter_2
  construct_circumcenter_2_object();

  Compute_squared_radius_2
  compute_squared_radius_2_object();

  /// @}
};

/**
 * \ingroup PkgTriangulation2TraitsClasses
 *
 * @class Robust_weighted_circumcenter_filtered_traits_2
 *
 * Upgrades the functors models of `Kernel::ConstructWeightedCircumcenter_2`,
 * and `Kernel::ComputeSquaredRadiusSmallestOrthogonalCircle_2`
 * that are provided by `K` with robust versions. This improved robustness
 * is achieved by using exact computations if the geometric configuration
 * is close to degenerate (e.g. almost coplanar points). The upgrade is completely
 * invisible from an interface point of view as the class `Robust_weighted_circumcenter_filtered_traits_2`
 * overrides the types and function objects provided by `K`.
 *
 * \tparam K must be a model of the `Kernel` concept.
 *
 * \sa `CGAL::Robust_circumcenter_filtered_traits_2<K>`
 * \sa \ref devman_robustness
 */
template<class K>
class Robust_weighted_circumcenter_filtered_traits_2
  : public Robust_circumcenter_filtered_traits_2<K>
{
public:
  /// \name Types
  /// @{

  /*!
    a model of `Kernel::ConstructWeightedCircumcenter_2`
  */
  typedef unspecified_type Construct_weighted_circumcenter_2;

  /*!
    a model of `Kernel::ComputeSquaredRadiusSmallestOrthogonalCircle_2`
  */
  typedef unspecified_type Compute_squared_radius_smallest_orthogonal_circle_2;

  /// @}

  /// \name Operations
  /// @{

  Construct_weighted_circumcenter_2
  construct_weighted_circumcenter_2_object();

  Compute_squared_radius_smallest_orthogonal_circle_2
  compute_squared_radius_smallest_orthogonal_circle_2_object();

  /// @}
};

}  // end namespace CGAL
