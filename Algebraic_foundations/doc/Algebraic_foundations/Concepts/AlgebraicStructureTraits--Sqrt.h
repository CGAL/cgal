
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction` providing the square root.

\cgalRefines{AdaptableUnaryFunction}

\sa `AlgebraicStructureTraits`

*/

class Sqrt {
public:

/// \name Types
/// @{

/*!
Is `AlgebraicStructureTraits::Type`.
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
returns \f$ \sqrt{x}\f$.
*/
result_type operator()(argument_type x) const;

/// @}

}; /* end Sqrt */

} /* end of namespace AlgebraicStructureTraits_ */
