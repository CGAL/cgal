namespace ArrDirectionalTraits {
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

\cgalRefines{AdaptableBinaryFunction}

\cgalHasModelsBegin
\cgalHasModels{ArrangementDirectionalXMonotoneTraits_2::Intersect_2}
\cgalHasModelsEnd

*/

class Intersect_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*! computes the intersections of `xc1` and `xc2` and inserts them <I>in an
 * ascending lexicographic \f$ xy\f$-order</I> into a range beginning at
 * `oi`. The type `OutputIterator` dereferences a `std::variant` of either the
 * type `pair<ArrDirectionalTraits::Point_2,
 * ArrDirectionalTraits::Multiplicity>` or the type
 * `ArrDirectionalTraits::X_monotone_curve_2`. An object of the former type
 * represents an intersection point with its multiplicity (in case the
 * multiplicity is undefined or unknown, it is set to \f$ 0\f$). An object of
 * the latter type representing an overlapping subcurve of `xc1` and `xc2`. The
 * overlapping subcurves are given the direction of `xc1` and `xc2` if their
 * directions are identical. Otherwise, the overlapping subcurves are given an
 * arbitrary direction. The operator returns a past-the-end iterator of the
 * destination range.
 */
OutputIterator operator()(ArrDirectionalTraits::X_monotone_curve_2 xc1,
                          ArrDirectionalTraits::X_monotone_curve_2 xc2,
                          Output_iterator& oi);

/// @}

}; /* end ArrDirectionalTraits::Intersect_2 */

}
