/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * `ArrangementHorizontalSideTraits_2` is an abstract concept. It generalizes all
 * concepts that handle curves that either reach or approach either the bottom
 * or top sizeds of the boundary of the parameter space.
 *
 * \cgalRefines `ArrangementBasicTraits_2`
 *
 * \cgalHasModel `CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, X, Y>`
 *
 * \sa `ArrangementVerticalSideTraits_2`
 */

class ArrangementHorizontalSideTraits_2 {
public:

  /// \name Categories
  /// @{
  /// @}

  /// \name Types
  /// @{
  /// @}

  /// \name Functor Types
  /// @{

  /// models the concept `ArrTraits::ParameterSpaceInX_2`.
  typedef unspecified_type Parameter_space_in_x_2;

  /// models the concept `ArrTraits::CompareXOnBoundary_2`.
  typedef unspecified_type Compare_x_on_boundary_2;

  /// models the concept `ArrTraits::CompareXNearBoundary_2`.
  typedef unspecified_type Compare_x_near_boundary_2;
  /// @}

  /// \name Accessing Functor Objects
  /// @{
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const;
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const;
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const;
  /// @}

}; /* end ArrangementHorizontalSideTraits_2 */
