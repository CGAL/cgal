/*!
 * \ingroup PkgArrangement2ConceptsTraits
 * \cgalConcept
 *
 * `ArrangementVerticalSideTraits_2` is an abstract concept. It generalizes all
 * concepts that handle curves that either reach or approach either the left
 * or right sizeds of the boundary of the parameter space.
 *
 * \cgalRefines `ArrangementBasicTraits_2`
 *
 * \cgalHasModel `CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, X, Y>`
 *
 * \sa `ArrangementVerticalSideTraits_2`
 */

class ArrangementVerticalSideTraits_2 {
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
  typedef unspecified_type Parameter_space_in_y_2;

  /// models the concept `ArrTraits::CompareXOnBoundary_2`.
  typedef unspecified_type Compare_y_on_boundary_2;

  /// models the concept `ArrTraits::CompareXNearBoundary_2`.
  typedef unspecified_type Compare_y_near_boundary_2;
  /// @}

  /// \name Accessing Functor Objects
  /// @{
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const;
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const;
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const;
  /// @}

}; /* end ArrangementHorizontalSideTraits_2 */
