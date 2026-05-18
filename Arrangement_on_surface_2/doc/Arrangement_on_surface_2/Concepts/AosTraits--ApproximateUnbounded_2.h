namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Approximate_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosApproximatePointTraits_2::Approximate_2}
 * \cgalHasModels{AosApproximateTraits_2::Approximate_2}
 * \cgalHasModels{AosApproximateUnboundedTraits_2::Approximate_2}
 * \cgalHasModelsEnd
 */
class ApproximateUnbounded_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! approximates a given \f$x\f$-monotone curve constrained to a bounding
   * box. It computes one or more sequences of approximate points that represent
   * the disconnected portions of a polyline that approximates `xcv` within the
   * bounding box `bbox`, and inserts them into output containers given through
   * the output iterator `oi`.  The first point of the first sequence and the
   * last point of the last sequence are always approximations of the endpoints
   * of the given curve.
   *
   * \param xcv The exact \f$x\f$-monotone curve.
   * \param error The error bound of the polyline approximation. This is the
   *        Hausdorff distance between the curve and the polyline that
   *        approximates the curve.
   * \param oi An output iterator for the output containers.
   * \param bbox the bounding box.
   * \param l2r A Boolean flag that indicates whether the curve direction is
   *        left to right.
   * \return The past-the-end iterator of the output container.
   *
   * \pre Dereferencing `oi` must yield an object the type of which is a
   *      container, where the value type of this container is
   *      `AosApproximateTraits_2::Approximate_point_2`.
   */
  template <typename OutputIterator>
  OutputIterator operator()(const X_monotone_curve_2& xcv, double error, OutputIterator oi,
                            const Bbox_2& bbox, bool l2r = true) const;

  /// @}
}; /* end AosTraits::Approximate_2 */

}
