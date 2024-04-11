/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * Several predicates are required to handle \f$ x\f$-monotone curves that
 * approach infinity and thus approach the boundary of the parameter
 * space. These predicates are sufficient to handle not only curves embedded in
 * an unbounded parameter space, but also curves embedded in a bounded parameter
 * space with open boundaries. Models of the concept
 * `ArrangementOpenBoundaryTraits_2` handle curves that approach the boundary of
 * a parameter space. This concept refines the concept
 * `ArrangementBasicTraits_2`. The arrangement template instantiated with a
 * traits class that models this concept can handle \f$ x\f$-monotone curves
 * that are unbounded in any direction. The concept
 * `ArrangementOpenBoundaryTraits_2`, nonetheless, also supports planar \f$
 * x\f$-monotone curves that reach the boundary of an open yet bounded parameter
 * space.
 *
 * An \f$ x\f$-monotone curve may be <I>closed</I>, in which case its endpoints
 * are representable as `Point_2` objects, or <I>open</I> at the boundary of the
 * parameter space. It can have one open end and one closed end (e.g., a
 * ray). The nature of the \f$ x\f$-monotone curves, whether they are expected
 * to be closed or not at any one of the four boundary-sides, is conveyed
 * through the definition of the four nested types `Left_side_category`,
 * `Right_side_category`, `Bottom_side_category`, and `Top_side_category`. If
 * some curves handled by a model of the concept
 * `ArrangementOpenBoundaryTraits_2` are expected to be open on the left, the
 * nested type `Left_side_category` must be convertible to
 * `CGAL::Arr_open_side_tag`. Similarly, if some curves handled by the concept
 * are expected to be open on the right, open at the bottom, or open at the top,
 * the corresponding nested type must be convertible to
 * `CGAL::Arr_open_side_tag`. A model of the concept
 * `ArrangementOpenBoundaryTraits_2` must have all the four categories
 * convertible to `CGAL::Arr_open_side_tag`.\cgalFootnote{We intend to introduce
 * more concepts that require only a subset of the categories to be convertible
 * to \cgalFootnoteCode{CGAL::Arr_open_side_tag}.} In this case the \dcel of the arrangement
 * instantiated with the model is initialized with an implicit bounding
 * rectangle. When the parameter space is bounded, it is the exact geometric
 * embedding of the implicit bounding rectangle.
 *
 * \cgalRefines{ArrangementBasicTraits_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arr_linear_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>}
 * \cgalHasModels{CGAL::Arr_algebraic_segment_traits_2<Coefficient>}
 * \cgalHasModels{CGAL::Arr_curve_data_traits_2<Tr,XData,Mrg,CData,Cnv>}
 * \cgalHasModels{CGAL::Arr_consolidated_curve_data_traits_2<Traits,Data>}
 * \cgalHasModelsEnd
 *
 * \sa `ArrangementBasicTraits_2`
 * \sa `ArrangementXMonotoneTraits_2`
 * \sa `ArrangementLandmarkTraits_2`
 * \sa `ArrangementTraits_2`
 */
class ArrangementOpenBoundaryTraits_2 {
public:

  /// \name Categories
  /// @{

  //! Must be convertible to `CGAL::Arr_open_side_tag`.
  typedef unspecified_type Left_side_category;

  //! Must be convertible to `CGAL::Arr_open_side_tag`.
  typedef unspecified_type Bottom_side_category;

  //! Must be convertible to `CGAL::Arr_open_side_tag`.
  typedef unspecified_type Top_side_category;

  //! Must be convertible to `CGAL::Arr_open_side_tag`.
  typedef unspecified_type Right_side_category;

  /// @}

  /// \name Functor Types
  /// @{

  /*! models the concept `ArrTraits::ParameterSpaceInX_2`.  Required only if the
   * traits class supports unbounded curves that approach the left or the right
   * sides (the `Left_side_category` or the `Right_side_category` categories are
   * convertible to `CGAL::Arr_open_side_tag`).
   */
  typedef unspecified_type Parameter_space_in_x_2;

  /*! models the concept `ArrTraits::CompareYNearBoundary_2`.
   * Required only if the traits class supports unbounded curves that approach
   * the left or the right sides (the `Left_side_category` or the
   * `Right_side_category` categories are convertible to
   * `CGAL::Arr_open_side_tag`).
   */
  typedef unspecified_type Compare_y_near_boundary_2;

  /*! models the concept `ArrTraits::ParameterSpaceInY_2`.
   * Required only if the traits class supports unbounded curves that approach
   * the bottom or the top sides (the `Bottom_side_category` or the
   * `Top_side_category` categories are convertible to
   * `CGAL::Arr_open_side_tag`).
   */
  typedef unspecified_type Parameter_space_in_y_2;

  /*! models the concept `ArrTraits::CompareXOnBoundaryOfCurveEnd_2`.  Required
   * only if the traits class supports unbounded curves that approach the bottom
   * or the top sides (the `Bottom_side_category` or the `Top_side_category`
   * categories are convertible to `CGAL::Arr_open_side_tag`).
   */
  typedef unspecified_type Compare_x_on_boundary_2;

  /*! models the concept `ArrTraits::CompareXNearBoundary_2`.  Required only if
   * the traits class supports unbounded curves that approach the bottom or the
   * top sides (the `Bottom_side_category` or the `Top_side_category` categories
   * are convertible to `CGAL::Arr_open_side_tag`).
   */
  typedef unspecified_type Compare_x_near_boundary_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{

  /*! */
  Parameter_space_in_x_2 parameter_space_in_x_2_object() const;

  /*! */
  Compare_y_near_boundary_2 compare_y_near_boundary_2_object() const;

  /*! */
  Parameter_space_in_y_2 parameter_space_in_y_2_object() const;

  /*! */
  Compare_x_on_boundary_2 compare_x_on_boundary_2_object() const;

  /*! */
  Compare_x_near_boundary_2 compare_x_near_boundary_2_object() const;

  /// @}

}; /* end ArrangementOpenBoundaryTraits_2 */
