namespace ArrTraits {
/*!
\ingroup PkgArrangementOnSurface2ConceptsFunctionObjects
\cgalConcept

\cgalRefines Functor

\cgalHasModel ArrangementTraits_2::Make_x_monotone_2

*/

class MakeXMonotone_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*! subdivides the input curve `c` into \f$ x\f$-monotone subcurves and
 * isolated points, and inserts the results into a range begining at the
 * given output iterator `oi`. The type `OutputIterator` dereferences a
 * `boost::variant` that wraps either an `ArrTraits::Point_2` object or an
 * `ArrTraits::X_monotone_curve_2` object. The operator returns a past-the-end
 * iterator for the output sequence.
 */
template <typename OutputIterator>
OutputIterator operator()(ArrTraits::Curve_2 c,
                          OutputIterator oi);

/// @}

}; /* end ArrTraits::MakeXMonotone_2 */

}
