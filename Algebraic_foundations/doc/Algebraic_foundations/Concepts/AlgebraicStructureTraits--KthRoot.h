
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableBinaryFunction` providing the k-th root.

\cgalRefines `AdaptableBinaryFunction`

\sa `FieldWithRootOf`
\sa `AlgebraicStructureTraits`

*/

class KthRoot {
public:

/// \name Types
/// @{

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type result_type;

/*!
Is int.
*/
typedef unspecified_type first_argument;

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type second_argument;

/// @}

/// \name Operations
/// @{

/*!
returns the \f$ k\f$-th root of \f$ x\f$.
\pre \f$ k \geq1\f$

*/
result_type operator()(int k, second_argument_type x);

/// @}

}; /* end KthRoot */

} /* end of namespace AlgebraicStructureTraits_ */
