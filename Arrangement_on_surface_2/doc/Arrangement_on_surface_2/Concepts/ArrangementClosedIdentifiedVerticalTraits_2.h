/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * `ArrangementClosedIdentifiedVerticalTraits_2` is an abstract concept. It
 * generalizes the concepts `ArrangementClosedLeftTraits_2`,
 * `ArrangementClosedRightTraits_2`, and
 * `ArrangementIdentifiedVerticalTraits_2`. (An "abstract" concept is a concept
 * that is useless on its own.) Only a combination of this concept and one or
 * more concepts that handle curves that either reach or approach the remaining
 * boundary sides (that is, borrom and top) are purposeful, and can have models.
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
