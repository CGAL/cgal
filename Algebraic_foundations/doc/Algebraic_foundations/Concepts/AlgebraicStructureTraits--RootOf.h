
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableFunctor` computes a real root of a square-free univariate
polynomial.

\cgalRefines `AdaptableFunctor`

\sa `FieldWithRootOf`
\sa `AlgebraicStructureTraits`

*/

class RootOf {
public:

/// \name Types
/// @{

/*!
Is `AlgebraicStructureTraits::Type`.
*/
typedef unspecified_type result_type;

/// @}

/// \name Operations
/// @{

/*!
returns the k-th real root of the univariate polynomial,
which is defined by the iterator range,
where begin refers to the constant term.
\pre The polynomial is square-free.
\pre The value type of the InputIterator is `AlgebraicStructureTraits::Type`.
*/
template<class InputIterator>
result_type operator() (int k, InputIterator begin, InputIterator end);

/// @}

}; /* end RootOf */

} /* end of namespace AlgebraicStructureTraits_ */
