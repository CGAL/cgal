namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Approximate_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosApproximateUnboundedTraits_2::Approximate_2}
 * \cgalHasModelsEnd
 */
class ApproximateUnbounded_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! approximates a given \f$x\f$-monotone curve clipped to a bounding box.
   * It computes a sequence of approximate points that represents a polyline
   * that approximates `xcv` clipped to the bounding box `bbox`, and inserts
   * the points into an output container given through the output iterator
   * `oi`, with one exception: every maximal portion of `xcv` that lies in the
   * vertical strip bounded by the left and right sides of `bbox`, but below
   * the bottom side (resp. above the top side) of `bbox`, is not removed.
   * Instead, it is represented by one or more points that lie in the strip
   * strictly below the bottom side (resp. above the top side) of `bbox`.
   * Whenever a point that lies inside `bbox` is adjacent in the sequence to a
   * point that lies outside `bbox`, the former lies on the boundary of `bbox`.
   * Portions of `xcv` that lie to the left or to the right of the strip are
   * removed.
   *
   * Consequently, (i) the maximal runs of consecutive points that lie inside
   * `bbox` (including its boundary) form polylines that approximate `xcv`
   * clipped to `bbox`, and (ii)
   * clamping the \f$y\f$-coordinates of all the points to the
   * \f$y\f$-range of `bbox` yields an approximation of the portion of `xcv`
   * that lies in the strip, where the parts below (resp. above) `bbox` are
   * replaced by the corresponding parts of the bottom (resp. top) side of
   * `bbox`. In particular, if `xcv` passes below or above `bbox`, the
   * sequence is not empty, although no point of `xcv` lies inside `bbox`. If
   * `xcv` does not intersect the strip, the sequence is empty.
   *
   * \param xcv The exact \f$x\f$-monotone curve.
   * \param error The error bound of the polyline approximation. This is the
   *        Hausdorff distance between the curve clipped to `bbox` and the
   *        polylines formed by the maximal runs of consecutive points that
   *        lie inside `bbox`.
   * \param oi An output iterator for the output container.
   * \param bbox the bounding box.
   * \param l2r A Boolean flag that indicates whether the curve direction is
   *        left to right.
   * \return The past-the-end iterator of the output container.
   *
   * \pre The bounding box `bbox` is not degenerate; that is, its width and
   *      height are positive.
   * \pre Dereferencing `oi` must yield an object the type of which is
   *      `AosApproximateTraits_2::Approximate_point_2`.
   */
  template <typename OutputIterator>
  OutputIterator operator()(const X_monotone_curve_2& xcv, double error, OutputIterator oi,
                            const Bbox_2& bbox, bool l2r = true) const;

  /// @}
}; /* end AosTraits::Approximate_2 */

}
