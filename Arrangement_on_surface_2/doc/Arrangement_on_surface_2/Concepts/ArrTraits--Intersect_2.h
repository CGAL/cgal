namespace ArrTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{Functor}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{ArrangementXMonotoneTraits_2::Intersect_2}
 * \cgalHasModelsEnd
 */
class Intersect_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! computes the intersections of two \f$x\f$-monotone curves and inserts the
   * result in ascending \f$xy\f$-lexicographic order into an output container
   * given through an output iterator. An intersection, if exists, is
   * represented by a discriminated union container that holds either an
   * intersection point along with its multiplicity or an overlapping
   * \f$x\f$-monotone subcurve. If the multiplicity is undefined or unknown, it
   * should be set to \f$0\f$.
   *
   * \param xc1 The first \f$x\f$-monotone curve.
   * \param xc2 The second \f$x\f$-monotone curve.
   * \param oi The output iterator that points at the output container.
   * \return The past-the-end iterator of the output container.
   *
   * A special case may occur when the left and right sides of the boundary of
   * the parameter space of the surface, the arrangement is embedded on, are
   * identified. An intersection point that lies on the identification curve,
   * between two \f$x\f$-monotone curves that intersect at their left and right
   * ends must be ignored.  Consider two \f$x\f$-monotone curves that intersect
   * at their left and right ends, respectively, at a point \f$p\f$ that lies on
   * the identification curve. If, for example, the number of intersections
   * between these two curves is greater than 1, the order of intersections is
   * nondeterministic.
   *
   * \pre Dereferencing `oi` must yield an object of type
   * `std::optional<std::variant<std::pair<%Point_2,ArrangementXMonotoneTraits_2::Multiplicity,X_monotone_curve_2>>`,
   * where `%Point_2` is a model of `ArrTraits::Point_2` and
   * `X_monotone_curve_2` is a model of `ArrTraits::XMonotoneCurve_2`.
   */
  template <typename OutputIterator>
  OutputIterator operator()(ArrTraits::X_monotone_curve_2 xc1,
                            ArrTraits::X_monotone_curve_2 xc2,
                            OutputIterator& oi);

/// @}

}; /* end ArrTraits::Intersect_2 */

}
