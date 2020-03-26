
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableFunctor` computes a square-free factorization
<I>up to a constant factor (utcf)</I> of a
`PolynomialTraits_d::Polynomial_d`.

A polynomial \f$ p\f$ is factored into square-free and pairwise coprime non-constant
factors \f$ g_i\f$ with multiplicities \f$ m_i\f$, such that
\f$ a \cdot p = g_1^{m_1} \cdot ... \cdot g_n^{m_n}\f$, where \f$ a\f$ is some constant factor.

The pairs \f$ (g_i,m_i)\f$ are written into the given output iterator.

The constant factor \f$ a\f$ is not computed.

This functor is well defined even though
`PolynomialTraits_d::Innermost_coefficient_type` may not be a
`UniqueFactorizationDomain`.

\cgalRefines Assignable
\cgalRefines CopyConstructible
\cgalRefines DefaultConstructible

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::SquareFreeFactorize`

*/

class PolynomialTraits_d::SquareFreeFactorizeUpToConstantFactor {
public:

/// \name Operations
/// @{

/*!
Computes the square-free factorization of \f$ p\f$ and returns the
past-the-end iterator of the written range.
\pre `std::iterator_traits< OutputIterator >::%value_type` must be constructible from `std::pair<PolynomialTraits_d::Polynomial_d,int>`.

*/
template<class OutputIterator>
OutputIterator operator()(PolynomialTraits_d::Polynomial_d p,
OutputIterator it);

/// @}

}; /* end PolynomialTraits_d::SquareFreeFactorizeUpToConstantFactor */

