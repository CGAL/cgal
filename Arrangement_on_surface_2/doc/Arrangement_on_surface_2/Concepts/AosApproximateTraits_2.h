/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * The concept `AosApproximateTraits_2` refines the concept
 * `AosApproximatePointTraits_2`. A model of this concept is able to
 * approximate a point and a curve (in addition to the ability to approximate the
 * coordinates of a point).
 *
 * \cgalRefines{AosApproximatePointTraits_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arr_circle_segment_traits_2<Kernel>}
 * \cgalHasModels{CGAL::Arr_conic_traits_2<RatKernel,AlgKernel,NtTraits>}
 * \cgalHasModels{CGAL::Arr_geodesic_arc_on_sphere_traits_2}
 * \cgalHasModels{CGAL::Arr_polyline_traits_2<SegmentTraits_2>}
 * \cgalHasModels{CGAL::Arr_segment_traits_2<Kernel>}
 * \cgalHasModelsEnd
 *
 * \sa `AosApproximatePointTraits_2`
 * \sa `draw()`
 */
class AosApproximateTraits_2 {
public:
  /// \name Types
  /// @{

  //! the approximate point.
  typedef unspecified_type Approximate_point_2;

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
