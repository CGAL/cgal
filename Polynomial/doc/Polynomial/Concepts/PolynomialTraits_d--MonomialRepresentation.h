
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `Functor` outputs the monomial representation of the given polynomial,
that is, it writes all non zero terms of the polynomial as
`std::pair<CGAL::Exponent_vector, PolynomialTraits_d::Innermost_coefficient_type>`
into the given output iterator.

\cgalRefines \ref Assignable
\cgalRefines \ref CopyConstructible
\cgalRefines \ref DefaultConstructible

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::ConstructPolynomial`

*/

class PolynomialTraits_d::MonomialRepresentation {
public:

/// \name Operations
/// @{

/*!
Writes the monom representation of \f$ p\f$ into the given output iterator \f$ it\f$.
\pre `std::iterator_traits< OutputIterator >::%value_type` must be constructible from `std::pair<CGAL::Exponent_vector, PolynomialTraits_d::Innermost_coefficient_type>`.

*/
template<class OutputIterator>
OutputIterator operator()(PolynomialTraits_d::Polynomial_d p,
OutputIterator it);

/// @}

}; /* end PolynomialTraits_d::MonomialRepresentation */

