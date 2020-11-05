
namespace RealEmbeddableTraits_ {

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction` computes the absolute value of a number.

\cgalRefines `AdaptableUnaryFunction`

\sa `RealEmbeddableTraits`

*/

class Abs {
public:

/// \name Types
/// @{

/*!
Is `RealEmbeddableTraits::Type`.
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
computes the absolute value of \f$ x\f$.
*/
result_type operator()(argument_type x);

/// @}

}; /* end Abs */

} /* end of namespace RealEmbeddableTraits_ */
