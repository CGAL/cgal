
namespace RealEmbeddableTraits_ {

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction`, returns true in case the argument is 0.

\cgalRefines{AdaptableUnaryFunction}

\sa `RealEmbeddableTraits`
\sa `AlgebraicStructureTraits_::IsZero`

*/

class IsZero {
public:

/// \name Types
/// @{

/*!
Type convertible to `bool`.
*/
typedef unspecified_type result_type;

/*!
Is `RealEmbeddableTraits::Type`.
*/
typedef unspecified_type argument_type;

/// @}

/// \name Operations
/// @{

/*!

returns true in case \f$ x\f$ is the zero element of the ring.
*/
result_type operator()(argument_type x);

/// @}

}; /* end IsZero */

} /* end of namespace RealEmbeddableTraits_ */
