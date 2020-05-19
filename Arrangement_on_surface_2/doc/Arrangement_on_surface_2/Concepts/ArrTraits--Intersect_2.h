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

/*! computes the intersections of `xc1` and `xc2` and inserts them <I>in an
 * ascending lexicographic \f$ xy\f$-order</I> into a range begining at
 * `oi`. The type `OutputIterator` dereferences a `boost::variant` of either the
 * type `pair<ArrTraits::Point_2,ArrTraits::Multiplicity>` or the type
 * `ArrTraits::X_monotone_curve_2`. An object of the former type represents an
 * intersection point with its multiplicity (in case the multiplicity is
 * undefined or unknown, it should be set to \f$ 0\f$). An object of the latter
 * type representing an overlapping subcurve of `xc1` and `xc2`. The operator
 * returns a past-the-end iterator of the destination range.
 */
OutputIterator operator()(ArrTraits::X_monotone_curve_2 xc1,
                          ArrTraits::X_monotone_curve_2 xc2,
                          Output_iterator& oi);

/// @}

}; /* end ArrTraits::Intersect_2 */

}
