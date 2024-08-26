namespace ArrTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Functor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{ArrangementTraits_2::Make_x_monotone_2}
 * \cgalHasModelsEnd
 */
class MakeXMonotone_2 {
public:
  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! subdivides an input curve into \f$x\f$-monotone subcurves and isolated
   * points, and inserts the results into an output container given through an
   * output iterator. An object in the output container is represented by a
   * discriminated union container that holds either a point or an
   * \f$x\f$-monotone curve.
   *
   * \param c The input curve.
   * \param oi The output iterator that points at the output container.
   * \return The past-the-end iterator of the output container.
   *
   * \pre Dereferencing `oi` must yield a polymorphic object of type
   * `std::variant<%Point_2, X_monotone_curve_2>`, where `%Point_2` is a model
   * of `ArrTraits::Point_2` and `X_monotone_curve_2` is a model of
   * `ArrTraits::XMonotoneCurve_2`.
   */
  template <typename OutputIterator>
  OutputIterator operator()(ArrTraits::Curve_2 c, OutputIterator oi);

  /// @}
}; /* end ArrTraits::MakeXMonotone_2 */

}
