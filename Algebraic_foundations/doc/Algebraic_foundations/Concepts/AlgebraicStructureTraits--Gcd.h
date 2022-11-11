
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableBinaryFunction` providing the gcd.

The greatest common divisor (\f$ gcd\f$) of ring elements \f$ x\f$ and \f$ y\f$ is the unique
ring element \f$ d\f$ (up to a unit) with the property that any common divisor of
\f$ x\f$ and \f$ y\f$ also divides \f$ d\f$. (In other words: \f$ d\f$ is the greatest lower bound
of \f$ x\f$ and \f$ y\f$ in the partial order of divisibility.) We demand the \f$ gcd\f$ to be
unit-normal (i.e.\ have unit part 1).

\f$ gcd(0,0)\f$ is defined as \f$ 0\f$, since \f$ 0\f$ is the greatest element with respect
to the partial order of divisibility. This is because an element \f$ a \in R\f$ is said to divide \f$ b \in R\f$, iff \f$ \exists r \in R\f$ such that \f$ a \cdot r = b\f$.
Thus, \f$ 0\f$ is divided by every element of the Ring, in particular by itself.

\cgalRefines `AdaptableBinaryFunction`

\sa `AlgebraicStructureTraits`

*/

class Gcd {
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
typedef unspecified_type first_argument;

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type second_argument;

/// @}

/// \name Operations
/// @{

/*!
returns \f$ gcd(x,y)\f$.
*/
result_type operator()(first_argument_type x,
second_argument_type y);

/*!
This operator is defined if `NT1` and `NT2` are `ExplicitInteroperable`
with coercion type `AlgebraicStructureTraits::Type`.
*/
template <class NT1, class NT2> result_type operator()(NT1 x, NT2 y);

/// @}

}; /* end Gcd */

} /* end of namespace AlgebraicStructureTraits_ */
