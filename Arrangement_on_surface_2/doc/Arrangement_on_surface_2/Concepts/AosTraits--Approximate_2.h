namespace AosTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{ApproximatePoint_2}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AosApproximatePointTraits_2::Approximate_2}
 * \cgalHasModels{AosApproximateTraits_2::Approximate_2}
 * \cgalHasModelsEnd
 */
class Approximate_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! obtains an approximation of `p`.
   */
  CGAL::Approximate_point_2 operator()(AosTraits::Point_2 p);

  /*! approximates a given \f$x\f$-monotone curve. It computes a sequence of
   * approximate points that represent an approximate polyline, and inserts
   * them into an output container given through an output iterator.  The
   * first and last points in the sequence are always approximations of the
   * endpoints of the given curve.
   *
   * \param xcv The exact \f$x\f$-monotone curve.
   * \param error The error bound of the polyline approximation. This is the
   *        Hausdorff distance between the curve and the polyline that
   *        approximates the curve.
   * \param oi An output iterator for the output container.
   * \param l2r A Boolean flag that indicates whether the curve direction is
   *        left to right.
   * \return The past-the-end iterator of the output container.
   *
   * \pre Dereferencing `oi` must yield an object of type
   *      `Arr_conic_traits_2::Approximate_point_2`.
   */
  template <typename OutputIterator>
  OutputIterator operator()(const X_monotone_curve_2& xcv, double error,
                            OutputIterator oi, bool l2r = true) const;

  /// @}
}; /* end AosTraits::Approximate_2 */

}
