/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * The concept `AosApproximatePointTraits_2` refines the basic traits concept
 * `AosBasicTraits_2`. A model of this concept is able to approximate a point.
 *
 * \cgalRefines{AosBasicTraits_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arr_circle_segment_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>}
 * \cgalHasModels{CGAL::Arr_geodesic_arc_on_sphere_traits_2}
 * \cgalHasModels{CGAL::Arr_linear_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_non_caching_segment_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_segment_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_polycurve_traits_2<GeometryTraits_2>}
 * \cgalHasModels{CGAL::Arr_polyline_traits_2<SegmentTraits_2>}
 * \cgalHasModels{CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>}
 * \cgalHasModelsEnd
 *
 * \sa `AosConstructXMonotoneCurveTraits_2`
 * \sa `AosXMonotoneTraits_2`
 * \sa `AosTraits_2`
 */
class AosApproximatePointTraits_2 {
public:
  /// \name Types
  /// @{

  //! the number type used to approximate point coordinates, e.g., double.
  typedef unspecified_type Approximate_number_type;

  /// @}

  /// \name Functor Types
  /// @{

  /// models the concept `AosTraits::Approximate_2`.
  typedef unspecified_type Approximate_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{

  ///
  Approximate_2 approximate_2_object() const;

  /// @}
}
