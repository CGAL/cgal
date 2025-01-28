
namespace RealEmbeddableTraits_ {

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction`, returns true in case the argument is negative.

\cgalRefines{AdaptableUnaryFunction}

\sa `RealEmbeddableTraits`

*/

class IsNegative {
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
returns true in case \f$ x\f$ is negative.
*/
result_type operator()(argument_type x);

/// @}

}; /* end IsNegative */

} /* end of namespace RealEmbeddableTraits_ */
