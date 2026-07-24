/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * The concept `AosConstructXMonotoneCurveTraits_2` refines the basic
 * traits concept `AosBasicTraits_2`. A model of this concept is able to
 * construct an \f$x\f$-monotone curve from two points.
 *
 * \cgalRefines{AosBasicTraits_2}
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
 * \sa `AosApproximatePointTraits_2`
 * \sa `AosXMonotoneTraits_2`
 * \sa `AosTraits_2`
 * \sa `AosConstructCurveTraits_2`
 */
class AosConstructXMonotoneCurveTraits_2 {
public:
  /// \name Functor Types
  /// @{

  /*! models the concept `AosTraits::ConstructXMonotoneCurve_2`.
   */
  typedef unspecified_type Construct_x_monotone_curve_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{

  ///
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2_object() const;

  /// @}
};
