
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableBinaryFunction`,
returns true if the first argument divides the second argument.

Integral division (a.k.a. exact division or division without remainder) maps
ring elements \f$ (n,d)\f$ to ring element \f$ c\f$ such that \f$ n = dc\f$ if such a \f$ c\f$
exists. In this case it is said that \f$ d\f$ divides \f$ n\f$.

This functor is required to provide two operators. The first operator takes two
arguments and returns true if the first argument divides the second argument.
The second operator returns \f$ c\f$ via the additional third argument.

\cgalRefines{AdaptableBinaryFunction}

\sa `AlgebraicStructureTraits`
\sa `AlgebraicStructureTraits_::IntegralDivision`

*/

class Divides {
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
typedef unspecified_type first_argument;

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type second_argument;

/// @}

/// \name Operations
/// @{

/*!
Computes whether \f$ d\f$ divides \f$ n\f$.
*/
result_type operator()(first_argument_type d,
second_argument_type n);

/*!

Computes whether \f$ d\f$ divides \f$ n\f$.
Moreover it computes \f$ c\f$ if \f$ d\f$ divides \f$ n\f$,
otherwise the value of \f$ c\f$ is undefined.

*/
result_type operator()(
first_argument_type d,
second_argument_type n,
AlgebraicStructureTraits::Type& c);

/// @}

}; /* end Divides */

} /* end of namespace AlgebraicStructureTraits_ */
