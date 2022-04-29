/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * The concept `ArrangementConstructCurveTraits_2` refines the basic
 * traits concept `ArrangementBasicTraits_2`. A model of this concept is able
 * to construct a curve from two points.
 *
 * \cgalRefines `ArrangementTraits_2`
 *
 * \cgalHasModel `CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>`
 * \cgalHasModel `CGAL::Arr_geodesic_arc_on_sphere_traits_2`
 * \cgalHasModel `CGAL::Arr_linear_traits_2<Kernel>`
 * \cgalHasModel `CGAL::Arr_non_caching_segment_traits_2<Kernel>`
 * \cgalHasModel `CGAL::Arr_segment_traits_2<Kernel>`
 * \cgalHasModel `CGAL::Arr_polyline_traits_2<SegmentTraits_2>`
 * \cgalHasModel `CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>`
 *
 * \sa `ArrangementConstructXMonotoneCurveTraits_2`, and
 *     `ArrangementTraits_2`
 */
class ArrangementConstructCurveTraits_2 {
public:
  /// \name Functor Types
  /// @{

  /*! models the concept `ArrTraits::ConstructCurve_2`.
   */
  typedef unspecified_type Construct_curve_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{

  /*!
   */
  Construct_curve_2 construct_curve_2_object() const;

  /// @}
};
