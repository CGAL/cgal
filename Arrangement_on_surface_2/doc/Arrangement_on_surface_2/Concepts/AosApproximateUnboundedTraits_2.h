/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * The concept `AosApproximateUnboundedTraits_2` refines the concept
 * `AosApproximateTraits_2`. A model of this concept is able to approximate a
 * curve constrained to a given bounding box (in addition to the ability to
 * approximate a point and a curve without constraints).
 *
 * \cgalRefines{AosApproximateTraits_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arr_linear_traits_2<Kernel>}
 * \cgalHasModelsEnd
 *
 * \sa `AosApproximateTraits_2`
 * \sa `draw()`
 */
class AosApproximateUnboundedTraits_2 {
public:
  /// \name Types
  /// @{

  /// @}

  /// \name Functor Types
  /// @{

  /// models the concept `AosTraits::ApproximateUnbounded_2`.
  typedef unspecified_type Approximate_2;

  /// @}

  /// \name Accessing Functor Objects
  /// @{

  ///
  Approximate_2 approximate_2_object() const;

  /// @}
}
