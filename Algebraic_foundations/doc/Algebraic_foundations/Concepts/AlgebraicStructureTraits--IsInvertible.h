
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction`,
returns `true` in case the argument is invertible in the ring.

\cgalRefines{AdaptableUnaryFunction}

\sa `AlgebraicStructureTraits`

*/

class IsInvertible {
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

returns true in case there exists an element \f$ y\f$ in the ring such that \f$ x * y \f$ is the one of the ring.
*/
result_type operator()(argument_type x);

/// @}

}; /* end IsInvetible */

} /* end of namespace AlgebraicStructureTraits_ */
