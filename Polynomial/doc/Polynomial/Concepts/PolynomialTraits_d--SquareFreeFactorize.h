
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `Functor` computes a square-free factorization
of a `PolynomialTraits_d::Polynomial_d`.

A polynomial \f$ p\f$ is factored into square-free and pairwise coprime non-constant
factors \f$ g_i\f$ with multiplicities \f$ m_i\f$ and a constant factor \f$ a\f$, such that
\f$ p = a \cdot g_1^{m_1} \cdot ... \cdot g_n^{m_n}\f$.

The pairs \f$ (g_i,m_i)\f$ are written into the given output iterator.

This functor is well defined if `PolynomialTraits_d::Polynomial_d` is a
`UniqueFactorizationDomain`.

\cgalRefines{Assignable,CopyConstructible,DefaultConstructible}

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::SquareFreeFactorizeUpToConstantFactor`
\sa `PolynomialTraits_d::MakeSquareFree`
\sa `PolynomialTraits_d::IsSquareFree`

*/

class PolynomialTraits_d::SquareFreeFactorize {
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
OutputIterator it,
PolynomialTraits_d::Innermost_coefficient_type& a);

/*!
As the first operator, just not computing the factor \f$ a\f$.
*/
template<class OutputIterator>
OutputIterator operator()(PolynomialTraits_d::Polynomial_d p,
OutputIterator it);

/// @}

}; /* end PolynomialTraits_d::SquareFreeFactorize */

