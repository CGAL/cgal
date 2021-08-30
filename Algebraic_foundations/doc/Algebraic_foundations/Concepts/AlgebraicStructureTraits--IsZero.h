
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction`, returns true in case the argument is the zero element of the ring.

\cgalRefines `AdaptableUnaryFunction`

\sa `AlgebraicStructureTraits`
\sa `RealEmbeddableTraits_::IsZero`

*/

class IsZero {
public:

/// \name Types
/// @{

/*!
Is `AlgebraicStructureTraits::Boolean`.
*/
typedef unspecified_type result_type;

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type argument_type;

/// @}

/// \name Operations
/// @{

/*!

returns true in case \f$ x\f$ is the zero element of the ring.
*/
result_type operator()(argument_type x) const;

/// @}

}; /* end IsZero */

} /* end of namespace AlgebraicStructureTraits_ */
