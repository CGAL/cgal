/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * The concept `ArrangementLandmarkTraits_2` refines the traits concepts
 * `ArrangementApproximateTraits_2` and
 * `ArrangementConstructXMonotoneCurveTraits_2`. The type of an arrangement
 * associated with the landmark point-location strategy (see
 * `CGAL::Arr_landmarks_point_location`) must be an instance of the
 * `CGAL::Arrangement_2<Traits,Dcel>` class template, where the Traits parameter
 * is substituted by a model of this concept.
 *
 * \cgalRefines{ArrangementApproximateTraits_2,ArrangementConstructXMonotoneCurveTraits_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>}
 * \cgalHasModels{CGAL::Arr_geodesic_arc_on_sphere_traits_2}
 * \cgalHasModels{CGAL::Arr_linear_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_non_caching_segment_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_segment_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_polycurve_traits_2<GeometryTraits_2>}
 * \cgalHasModels{CGAL::Arr_polyline_traits_2<SegmentTraits_2>}
 * \cgalHasModels{CGAL::Arr_rational_function_traits_2}
 * \cgalHasModelsEnd
 *
 * \sa `ArrangementXMonotoneTraits_2` and
 * `ArrangementTraits_2<AlgebraicKernel_d_1>`
 */
class ArrangementLandmarkTraits_2 {
public:

  /// \name Types
  /// @{
  /// @}

  /// \name Functor Types
  /// @{
  /// @}

  /// \name Accessing Functor Objects
  /// @{
  /// @}

};
