/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * The concept `ArrangementConstructXMonotoneCurveTraits_2` refines the basic
 * traits concept `ArrangementBasicTraits_2`. A model of this concept is able to
 * construct an \f$ x\f$-monotone curve from two points.
 *
 * \cgalRefines{ArrangementBasicTraits_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>}
 * \cgalHasModels{CGAL::Arr_geodesic_arc_on_sphere_traits_2}
 * \cgalHasModels{CGAL::Arr_linear_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_non_caching_segment_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_segment_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_polyline_traits_2<SegmentTraits_2>}
 * \cgalHasModels{CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>}
 * \cgalHasModelsEnd
 *
 * \sa `ArrangementApproximateTraits_2`,
 *     `ArrangementXMonotoneTraits_2`
 *     `ArrangementTraits_2`, and
 *     `ArrangementConstructCurveTraits_2`.
 */
class ArrangementConstructXMonotoneCurveTraits_2 {
public:
  /// \name Functor Types
  /// @{

  /*! models the concept `ArrTraits::ConstructXMonotoneCurve_2`.
   */
  typedef unspecified_type Construct_x_monotone_curve_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{

  /*!
   */
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const;

  /// @}
};
