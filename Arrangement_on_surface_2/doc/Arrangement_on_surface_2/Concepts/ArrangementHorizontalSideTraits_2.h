/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * `ArrangementHorizontalSideTraits_2` is an abstract concept. It generalizes
 * all concepts that handle curves that either reach or approach either the
 * bottom or top sizeds of the boundary of the parameter space. (An "abstract"
 * concept is a concept that is useless on its own.) Only a combination of this
 * concept and one or more concepts that handle curves that either reach or
 * approach the remaining boundary sides (that is, left and right) are
 * purposeful, and can have models.
 *
 * \cgalRefines{ArrangementBasicTraits_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arr_linear_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>}
 * \cgalHasModels{CGAL::Arr_algebraic_segment_traits_2<Coefficient>}
 * \cgalHasModels{CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, X, Y>}
 * \cgalHasModelsEnd
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

  /// models the concept `ArrTraits::ParameterSpaceInY_2`.
  typedef unspecified_type Parameter_space_in_y_2;

  /// models the concept `ArrTraits::CompareXOnBoundaryOfCurveEnd_2`.
  typedef unspecified_type Compare_x_on_boundary_2;

  /// models the concept `ArrTraits::CompareXNearBoundary_2`.
  typedef unspecified_type Compare_x_near_boundary_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const;
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const;
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const;
  /// @}

}; /* end ArrangementHorizontalSideTraits_2 */
