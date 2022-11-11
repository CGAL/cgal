namespace ArrTraits {

/*! \ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines Functor
 *
 * \cgalHasModel ArrangementXMonotoneTraits_2::Intersect_2
 */
class Intersect_2 {
public:

  /// \name Operations
  /// A model of this concept must provide:
  /// @{

  /*! computes the intersections of `xc1` and `xc2` and writes them <I>in an
   * ascending lexicographic \f$xy\f$-order</I> into a range beginning at
   * `oi`. The type `OutputIterator` must dereference a polymorphic object of
   * type `boost::variant` that wraps objects of type either type
   * `pair<ArrTraits::Point_2, ArrTraits::Multiplicity>` or
   * `ArrTraits::X_monotone_curve_2`. An object of the former type represents an
   * intersection point with its multiplicity (in case the multiplicity is
   * undefined or unknown, it should be set to \f$0\f$). An object of the latter
   * type represents an overlapping subcurve of `xc1` and `xc2`. The operator
   * returns a past-the-end iterator of the destination range.
   *
   * A special case may occur when the parameter space of the surface, the
   * arrangement is embedded on, is identified on the left and right sides of
   * the boundary. An intersection point that lies on the identification curve,
   * between two \f$x\f$-monotone curves that intersect at their left and right
   * ends must be ignored.  Consider two \f$x\f$-monotone curves that intersect
   * at their left and right ends, respectively, at a point \f$p\f$ that lies on
   * the identification curve. If, for example, the number of intersections
   * between these two curves is greater than 1, the order of intersections is
   * non-deterministic.
   */
  template <typename OutputIterator>
  OutputIterator operator()(ArrTraits::X_monotone_curve_2 xc1,
                             ArrTraits::X_monotone_curve_2 xc2,
                             OutputIterator& oi);

/// @}

}; /* end ArrTraits::Intersect_2 */

}
