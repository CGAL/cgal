namespace ArrTraits {
/*!
\ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
\cgalConcept

\cgalRefines Functor

\cgalHasModel ArrangementXMonotoneTraits_2::Intersect_2

*/

class Intersect_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*! computes the intersections of `xc1` and `xc2` and writes them <I>in an
 * ascending lexicographic \f$ xy\f$-order</I> into the output iterator
 * `oi`. The type of a value written into `oi` must be convertible to
 * `CGAL::Object`. The value itself wraps either a value of type
 * `pair<ArrTraits::Point_2,ArrTraits::Multiplicity>` or a value of type
 * `ArrTraits::X_monotone_curve_2`. A value of of the former type represents an
 * intersection point with its multiplicity; in case the multiplicity is
 * undefined or unknown, it should be set to \f$ 0\f$). A value of the latter
 * type representing an overlapping subcurve of `xc1` and `xc2`. The operator
 * returns a past-the-end iterator for the output sequence.
 *
 * A special case may occur when the parameter space of the surface, the
 * arrangement is embedded on, is identified on the left and right sides of the
 * boundary. An intersection point that lies on the identification curve,
 * between two \f$X\f$-monotone curves that intersect at their left and right
 * ends must be ignored.  Consider two \f$X\f$-monotone curves that intersect at
 * their left and right ends, respectively, at a point \f$p\f$ that lies on the
 * identification curve. If, for example, the number of intersections between
 * these two curves is greater than 1, the order of intersections is
 * non-deterministic.
 */
Output_iterator operator()(ArrTraits::X_monotone_curve_2 xc1,
ArrTraits::X_monotone_curve_2 xc2,
Output_iterator& oi);

/// @}

}; /* end ArrTraits::Intersect_2 */

}
