/*! \ingroup PkgArrangementOnSurface2ConceptsTraits
 * \cgalConcept
 *
 * The concept `AosApproximateUnboundedTraits_2` refines the concept
 * `AosApproximateTraits_2`. A model of this concept is able to approximate an
 * \f$x\f$-monotone curve clipped to a given bounding box, where the portions
 * of the curve that lie below or above the box are represented by points that
 * lie below or above the box, respectively (in addition to the ability to
 * approximate a point and a curve without constraints).
 *
 * \cgalRefines{AosApproximateTraits_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arr_linear_traits_2<Kernel>}
 * \cgalHasModelsEnd
 *
 * \sa `AosApproximateTraits_2`
 * \sa \link PkgArrangementOnSurface2Draw `CGAL::draw()`\endlink
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
