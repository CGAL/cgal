
namespace RealEmbeddableTraits_ {

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction`, returns true in case the argument is positive.

\cgalRefines{AdaptableUnaryFunction}

\sa `RealEmbeddableTraits`

*/

class IsPositive {
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
returns true in case \f$ x\f$ is positive.
*/
result_type operator()(argument_type x);

/// @}

}; /* end IsPositive */

} /* end of namespace RealEmbeddableTraits_ */
