/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * `ArrangementClosedIdentifiedVerticalTraits_2` is an abstract concept. It
 * generalizes the concepts `ArrangementClosedLeftTraits_2`,
 * `ArrangementClosedRightTraits_2`, `ArrangementIdentifiedLeftTraits_2` and
 * `ArrangementIdentifiedRightTraits_2`.
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

  /// models the concept `ArrTraits::CompareXOnBoundary_2`.
  typedef unspecified_type Compare_y_on_boundary_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{
  Compare_y_on_boundary_2 compare_y_on_boundary_2_object() const;
  /// @}

}; /* end ArrangementHorizontalSideTraits_2 */
